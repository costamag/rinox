/* rinox: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file array_area_resyn.hpp
  \brief Resynthesis using multiple output gates

  \author Andrea Costamagna
*/

#ifndef RINOX_SYNTHESIS_LUT_SPIN_DEC_SOLO_HPP_
#define RINOX_SYNTHESIS_LUT_SPIN_DEC_SOLO_HPP_

#include <fmt/format.h>

#include <chrono>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <optional>

#include "utils.hpp"
#include "../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../windowing/window_manager.hpp"
#include "../../windowing/window_simulator.hpp"
#include "../../boolean/support_minimizer.hpp"
#include "../../boolean/spfd.hpp"

namespace rinox {
namespace synthesis {
namespace lut {

template<typename Information, uint32_t MaxNumSpins>
class SpinSystem
{
  struct Spin{
    uint32_t name;
    Information info;
    Information care;
    int value;
  };

public:
  void AddSpin(Information const& info, Information const& care, uint32_t name)
  {
    int n = spins_.size();
    for ( int i = 0; i < n; ++i )
    {
      int Jp = kitty::count_ones( info ^ spins_[i].info );
      int Jn = kitty::count_ones( ~info ^ spins_[i].info );
      J[n][i] = J[i][n] = (Jp-Jn)/2;
    }

    Spin s{name, info, care, 0};
    spins_.push_back(s);
  }

  void clear()
  {
    spins_.clear();
  }

  kitty::ternary_truth_table<Information> GetFunction()
  {
    std::vector<int> best_state(spins_.size(), -1);
    int best_energy = std::numeric_limits<int>::max();

    int n = spins_.size();
    for (int state = 0; state < (1u << n); ++state)
    {
      for (int i = 0; i < n; ++i)
        spins_[i].value = ( ( state >> i ) & 0x1 ) ? 1 : -1;
      int energy = ComputeEnergy();
      if (energy < best_energy)
      {
        best_energy=energy;
        for (int i = 0; i < n; ++i)
          best_state[i] = spins_[i].value;
      }
    }
    Information care;
    Information bits;
    for (int i = 0; i < n; ++i)
      care |= spins_[i].care;
    for (int i = 0; i < n; ++i)
      bits |= best_state[i] > 0 ? spins_[i].care&spins_[i].info : spins_[i].care&(~spins_[i].info);

    kitty::ternary_truth_table<Information> func(bits, care);
    return func;
  }

  int ComputeEnergy() const
  {
    int res = 0;
    const int n = static_cast<int>(spins_.size());
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        res -= J[i][j] * spins_[i].value * spins_[j].value;
    return res;
  }

private:
  int J [MaxNumSpins][MaxNumSpins];
  math::static_vector<Spin, MaxNumSpins> spins_;
};

template <class Ntk, typename Manager, typename Simulator,
          typename Params = DefaultLutSynthParams>
class SpinDecSolo {
  static constexpr uint32_t kMaxNumOutputs = Ntk::max_num_outputs;
  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  // the functional information is provided in a larger truth table
  using Signature = kitty::static_truth_table<Params::kNumVarsCube>;
  // the decompsition happens in a smaller Boolean space
  using Function = kitty::static_truth_table<Params::kMaxCutsSize>;
  using Arrivals = analyzers::trackers::arrival_times_tracker<Ntk>;
  using Reminders = SubProblems<Params>;
  using Cut = math::static_vector<Signal, Params::kMaxCutsSize>;
  using List = synthesis::lut::LutList<Params::kMaxInpsSolo, kMaxNumOutputs, Params::kMaxCutsSize>;
  using SuppMinimizer = rinox::boolean::support_minimizer<Params::kMaxCutsSize, Params::kExactSuppMin>;
  using Spfds = boolean::spfds_manager<Function, 1u << Params::kMaxInpsSolo, kMaxNumOutputs>;


 public:
  SpinDecSolo(Ntk&, Manager&, Simulator&, LutSynthStats&);
  std::optional<Signal> Run( Cut const&, Signature const&, Signature const&);
  std::optional<typename List::Signal> RunRecursion(kitty::ternary_truth_table<Function> const&, std::vector<typename List::Signal> const&);
  math::static_vector<typename List::Signal, Params::kMaxInpsSolo> ChooseVars(std::vector<typename List::Signal>const&, kitty::ternary_truth_table<Function> const&,size_t);
  kitty::ternary_truth_table<Function> ExtractOneSimple();
  kitty::ternary_truth_table<Function> ExtractOneSpins(std::vector<typename List::Signal>const&);
  Signal Insert(std::vector<Signal> &, typename List::LutNode & n, typename List::Index & index);
  Function Generalize(Function const&, Function const&, std::vector<typename List::Signal> const&);

 private:
  Ntk& ntk_;
  List list_;
  Manager& manager_;
  Simulator& simulator_;
  LutSynthStats& st_;
  Reminders reminders_;
  std::array<Function, Params::kMaxCutsSize> proj_fns_;
  SuppMinimizer supp_minimizer_;
  Spfds spfds_;
  SpinSystem<Function, 1u << Params::kMaxInpsSolo> spins_;
};

template <class Ntk, typename Manager, typename Simulator, typename Params>
SpinDecSolo<Ntk, Manager, Simulator, Params>::SpinDecSolo(Ntk& ntk, Manager& manager, Simulator& simulator, LutSynthStats& st)
  : ntk_(ntk), manager_(manager), simulator_(simulator), st_(st) {
  for (int i = 0; i < Params::kMaxCutsSize; ++i) {
    kitty::create_nth_var(proj_fns_[i], i);
  }
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
auto SpinDecSolo<Ntk, Manager, Simulator, Params>::Run( Cut const& cut, Signature const& onset, Signature const& careset )->std::optional<Signal>
{
  std::vector<Signal> signals;
  reminders_.clear();
  std::vector<Signature const *> tts;
  std::vector<typename List::Signal> supp;

  list_.Reset(Params::kMaxNumGates);
  for ( uint32_t i = 0; i < cut.size(); ++i )
  {
    tts.push_back( &simulator_.Get(cut[i]) );
    supp.push_back(list_.MakeSignal(list_.CreateInput(&proj_fns_[i]),0));
    signals.push_back(cut[i]);
  }

  auto func = rinox::synthesis::lut::Interpolate<Params::kMaxCutsSize, Params::kNumVarsCube>( tts, onset, careset );
  auto res = RunRecursion(func, supp);
  if (!res) return std::nullopt;
  
  std::optional<Signal> fnew;
  list_.ForeachGate([&](auto & n, auto & index){
    fnew = std::make_optional(Insert( signals, n, index ));
  });
  return fnew;
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
auto SpinDecSolo<Ntk, Manager, Simulator, Params>::RunRecursion(kitty::ternary_truth_table<Function> const& func, std::vector<typename List::Signal> const& supp)-> std::optional<typename List::Signal>
{
  // minimize the support
  auto reminder = func;
  auto support = supp;
  std::vector<uint32_t> indices;
  for (auto const x : support)
    indices.push_back( list_.GetIndex(x) );

  supp_minimizer_.run(reminder, indices, support);

  // if smaller than the fanin size synthesize it
  if ( support.size() <= Params::kMaxInpsSolo )
  {
    math::static_vector<typename List::Signal, Params::kMaxInpsSolo> fanin;
    for (auto x : support) fanin.push_back(x);
    math::static_vector<kitty::ternary_truth_table<Function>, kMaxNumOutputs> funcs;
    funcs.push_back(func);
    return list_.CreateInternal(fanin, funcs);
  }
  // else recursive step
  auto supp1 = ChooseVars( support, func, Params::kMaxInpsSolo-1 );
  auto reminder1 = ExtractOneSpins(support);
  if ( kitty::count_ones(reminder1._care) < kitty::count_ones(func._care) )
  {
    auto f1 = RunRecursion( reminder1, support );
    if (f1)
    {
      supp1.push_back(*f1);
      math::static_vector<kitty::ternary_truth_table<Function>, kMaxNumOutputs> funcs;
      funcs.push_back(func);
      return list_.CreateInternal(supp1, funcs);
    }
  }
  return std::nullopt;
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
auto SpinDecSolo<Ntk, Manager, Simulator, Params>::ChooseVars(std::vector<typename List::Signal>const& supp0, kitty::ternary_truth_table<Function> const& func, size_t num_vars) -> math::static_vector<typename List::Signal, Params::kMaxInpsSolo>
{
  spfds_.init( {&func._bits}, {&func._care} );
  math::static_vector<typename List::Signal, Params::kMaxInpsSolo> res;
  uint32_t best_cost = spfds_.get_num_edges();
  while (!spfds_.is_covered() && (res.size() < num_vars))
  {
    std::optional<typename List::Signal> best;
    for (auto const& s : supp0)
    {
      auto const cost = spfds_.evaluate( list_.GetSimBits(s) );
      if (cost < best_cost)
      {
        best_cost = cost;
        best = std::make_optional(s);
      }
    }
    if (best)
    {
      res.push_back(*best);
      spfds_.update(list_.GetSimBits(*best));
    }
    else
      break;
  }
  return res;
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
auto SpinDecSolo<Ntk, Manager, Simulator, Params>::ExtractOneSimple() -> kitty::ternary_truth_table<Function>
{
  kitty::ternary_truth_table<Function> func;
  auto const& spfd = spfds_.Get(0);
  for ( uint32_t m = 0; m < spfd.num_masks_; ++m )
  {
    if (!spfd.is_killed(m))
    {
      auto const num_1s = kitty::count_ones(spfd.func_[1]&spfd.masks_[m]);
      auto const num_0s = kitty::count_ones(spfd.func_[0]&spfd.masks_[m]);
      func._bits |= ( num_1s >= num_0s ) ? spfd.func_[0]&spfd.masks_[m] : spfd.func_[1]&spfd.masks_[m];
      func._care |= spfd.masks_[m];
    }
  }
  return func;
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
auto SpinDecSolo<Ntk, Manager, Simulator, Params>::ExtractOneSpins(std::vector<typename List::Signal> const& support) -> kitty::ternary_truth_table<Function>
{
  spins_.clear();
  auto const& spfd = spfds_.Get(0);
  for ( uint32_t m = 0; m < spfd.num_masks_; ++m )
  {
    if (!spfd.is_killed(m))
    {
      Function const func = Generalize( spfd.func_[1], spfd.masks_[m], support );
      spins_.AddSpin( func, spfd.masks_[m], m );
    }
  }

  kitty::ternary_truth_table<Function> res = spins_.GetFunction();
  return res;
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
auto SpinDecSolo<Ntk, Manager, Simulator, Params>::Insert(std::vector<Signal> & signals, typename List::LutNode & n, typename List::Index & index ) -> Signal
{
  if (auto internal = std::get_if<typename List::Internal>(&n)) {
    // internal is a pointer to the Internal object stored in the variant
    auto& fanins = internal->fanins_;
    // use fanins as needed
    std::vector<Function const *> inps;
    std::vector<Signal> supp;
    for ( auto s : internal->fanins_ )
    {
      inps.push_back( &list_.GetSimBits( s ) );
      supp.push_back( signals[list_.GetIndex( s )] );
    }
    auto tts = Interpolate<Function, kMaxNumOutputs>( inps, internal->sims_ );
    auto fnew = ntk_.create_node( supp, tts );
    auto sims = ntk_.compute( ntk_.get_node( fnew ), inps );
    for ( int i = 0; i < sims.size(); ++i )
    {
      internal->sims_[i]._bits=sims[i];
      internal->sims_[i]._care|=~internal->sims_[i]._care;
      signals.push_back( Signal(ntk_.get_node(fnew), i) );
    }
    return signals.back();
  } else {
    throw std::runtime_error("Input not internal");
  }
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
auto SpinDecSolo<Ntk, Manager, Simulator, Params>::Generalize(Function const& bits, Function const& care, std::vector<typename List::Signal> const& support) -> Function
{
  struct TreeNode{
    typename List::Signal X;
    uint32_t T;
    uint32_t E;
  };

  math::static_vector<TreeNode, (1u << Params::kMaxInpsSolo)+2u> bdd;

  bdd.push_back({0,0,0});
  bdd.push_back({1,1,1});

  std::function<uint32_t(Function const&, Function const&)> synthesize = [&]( Function const& onset, Function const& careset ){
    kitty::ternary_truth_table<Function> fun(onset, careset);
    // chose the variable
    std::optional<typename List::Signal> best;

    auto const Ny = kitty::count_ones(careset & onset);
    int best_info = Ny*( kitty::count_ones(careset) - Ny );
    for (auto const s : support)
    {
      auto const& ttX = list_.GetSimBits( s );
      auto const care0 = careset & ~ttX;
      auto const care1 = careset & ttX;
      auto const Nx0 = kitty::count_ones( care0 & bits );
      auto const Nx1 = kitty::count_ones( care1 & bits );

      int info0 = Nx0*( kitty::count_ones(care0) - Nx0 );
      int info1 = Nx1*( kitty::count_ones(care1) - Nx1 );
      int info = info0 + info1;

      if (info<=best_info)
      {
        best_info = info;
        best = s;
      }
    }
    auto const& ttX = list_.GetSimBits( *best );
    auto const& ttT = careset & ttX;
    auto const& ttE = careset &~ttX;
    uint32_t T;
    uint32_t E;
    if ( kitty::count_ones( ttT & onset ) == 0 )
      T = 0u;
    else if ( kitty::equal( ttT & onset, ttT ) )
      T = 1u;
    else
      T = synthesize( onset, ttT );

    if ( kitty::count_ones( ttE & onset ) == 0 )
      E = 0u;
    else if ( kitty::equal( ttE & onset, ttE ) )
      E = 1u;
    else
      E = synthesize( onset, ttE );
    
    uint32_t index = bdd.size();
    TreeNode n { *best, T, E };
    bdd.push_back( n );
    return index;
  };

  synthesize(bits, care);

  math::static_vector<Function, (1u << Params::kMaxInpsSolo)+2u> sims;
  sims.push_back(bits^bits);
  sims.push_back(~sims[0]);
  Function res;
  for ( int i = 2; i < bdd.size(); ++i )
  {
    auto const& n = bdd[i];
    auto const& ttX = list_.GetSimBits( n.X );
    auto const& ttT = sims[n.T];
    auto const& ttE = sims[n.E];
    sims.push_back( (ttX & ttT) | (~ttX & ttE) );
  }
  return sims[sims.size()-1];
}


}  // namespace lut
}  // namespace synthesis
}  // namespace rinox

#endif