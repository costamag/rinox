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
#include "../../math/ising.hpp"
#include "../../synthesis/partial_specs.hpp"

namespace rinox {
namespace synthesis {
namespace lut {

template <class Ntk, typename Manager, typename Simulator,
          typename Params = DefaultLutSynthParams>
class SpinDecSolo {
  static constexpr uint32_t kMaxNumOutputs = Ntk::kMaxNumOutputs;
  static constexpr uint32_t kMaxNumSpins = 1u << Params::kMaxInpsSolo;

  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  // the functional information is provided in a larger truth table
  using Signature = kitty::static_truth_table<Params::kNumVarsCube>;
  // the decompsition happens in a smaller Boolean space
  using Function = kitty::static_truth_table<Params::kMaxCutsSize>;
  using Arrivals = analyzers::trackers::arrival_times_tracker<Ntk>;
  using Reminders = math::static_vector<kitty::ternary_truth_table<Function>, kMaxNumOutputs>;
  using Cut = math::static_vector<Signal, Params::kMaxCutsSize>;
  using List = synthesis::lut::LutList<Params::kMaxInpsSolo, kMaxNumOutputs, Params::kMaxCutsSize>;
  using SuppMinimizer = rinox::boolean::support_minimizer<Params::kMaxCutsSize, Params::kExactSuppMin>;
  using Spfds = boolean::spfds_manager<Function, 1u << Params::kMaxInpsSolo, kMaxNumOutputs>;
  using Learner = synthesis::partial::Learner<List, Function, Params::kMaxInpsSolo>;
  using ListNodeSignals = math::static_vector<typename List::Signal, kMaxNumOutputs>;
  using Spins = math::SpinSystem<Function, kMaxNumSpins, kMaxNumOutputs>;

  struct RemindersWrap {
    using type = std::conditional_t<kMaxNumOutputs==1, 
                                    std::array<Reminders, 2u>,
                                    Reminders>;
  };
  using RemindersType = typename RemindersWrap::type;


 public:
  SpinDecSolo(Ntk&, Manager&, Simulator&, LutSynthStats&);
  std::optional<Signal> Run( Cut const&, math::static_vector<Signature, kMaxNumOutputs> const&, math::static_vector<Signature, kMaxNumOutputs> const&);
  std::optional<ListNodeSignals> RunRecursion( Reminders const&, std::vector<typename List::Signal> const&);
  math::static_vector<typename List::Signal, Params::kMaxInpsSolo> ChooseVars(std::vector<typename List::Signal>const&, Reminders const&, size_t);
  RemindersType GetFunctions();
  std::optional<RemindersType> ExtractReminders(std::vector<typename List::Signal>const&);
  Signal Insert(std::vector<Signal> &, typename List::LutNode & n, typename List::Index & index);

 private:
  Ntk& ntk_;
  List list_;
  Manager& manager_;
  Simulator& simulator_;
  LutSynthStats& st_;
  std::array<Function, Params::kMaxCutsSize> proj_fns_;
  SuppMinimizer supp_minimizer_;
  Spfds spfds_;
  Spins spins_;
  Learner learner_;
};

template <class Ntk, typename Manager, typename Simulator, typename Params>
SpinDecSolo<Ntk, Manager, Simulator, Params>::SpinDecSolo(Ntk& ntk, Manager& manager, Simulator& simulator, LutSynthStats& st)
  : ntk_(ntk), manager_(manager), simulator_(simulator), st_(st), learner_(list_) {
  for (int i = 0; i < Params::kMaxCutsSize; ++i) {
    kitty::create_nth_var(proj_fns_[i], i);
  }
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
auto SpinDecSolo<Ntk, Manager, Simulator, Params>::Run( Cut const& cut, math::static_vector<Signature, kMaxNumOutputs> const& onset, math::static_vector<Signature, kMaxNumOutputs> const& careset )->std::optional<Signal>
{
  std::vector<Signal> signals;
  std::vector<Signature const *> tts;
  std::vector<typename List::Signal> supp;

  list_.Reset(Params::kMaxNumGates);
  for ( uint32_t i = 0; i < cut.size(); ++i )
  {
    tts.push_back( &simulator_.Get(cut[i]) );
    supp.push_back(list_.MakeSignal(list_.CreateInput(&proj_fns_[i]),0));
    signals.push_back(cut[i]);
  }

  auto const funcs = rinox::synthesis::lut::Interpolate<Params::kMaxCutsSize, Params::kNumVarsCube, kMaxNumOutputs>( tts, onset, careset );

  auto res = RunRecursion(funcs, supp);
  if (!res) return std::nullopt;
  
  std::optional<Signal> fnew;
  list_.ForeachGate([&](auto & n, auto & index){
    fnew = std::make_optional(Insert( signals, n, index ));
  });
  return fnew;
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
auto SpinDecSolo<Ntk, Manager, Simulator, Params>::RunRecursion(Reminders const& funcs, std::vector<typename List::Signal> const& supp)-> std::optional<ListNodeSignals>
{
  // minimize the support
  auto targets = funcs;
  auto support = supp;
  std::vector<uint32_t> indices;
  for (auto const x : support)
    indices.push_back( list_.GetIndex(x) );

  // TODO: implement the support minimizer for multiple outputs
  supp_minimizer_.template Run<kMaxNumOutputs>(targets, indices, support);

  bool is_done = ( funcs.size() == 1u ) && ( support.size() <= Params::kMaxInpsSolo );
  if constexpr ( kMaxNumOutputs == 2 )
    is_done |= ( ( funcs.size() == 2u ) && ( support.size() <= Params::kMaxInpsDual ) );
  if ( is_done )
  {
    math::static_vector<typename List::Signal, Params::kMaxInpsSolo> fanin;
    for (auto x : support) fanin.push_back(x);
    return list_.CreateInternal(fanin, funcs);
  }

  std::optional<math::static_vector<typename List::Signal, kMaxNumOutputs>> res;
  auto supp_top = ChooseVars( support, targets, Params::kMaxInpsSolo-2 );
  auto reminders = ExtractReminders(support);
  if (!reminders) return std::nullopt;
  if constexpr (kMaxNumOutputs == 1)
  {
    auto f1 = RunRecursion( (*reminders)[0], support );
    auto f2 = RunRecursion( (*reminders)[1], support );
    if (f2 && f1)
    {
      supp_top.push_back((*f1)[0]);
      supp_top.push_back((*f2)[0]);
      return list_.CreateInternal(supp_top, funcs);
    }
  }
  else if constexpr (kMaxNumOutputs == 2)
  {
    auto fs = RunRecursion( *reminders, support );
    if (fs)
    {
      supp_top.push_back((*fs)[0]);
      supp_top.push_back((*fs)[1]);
      return list_.CreateInternal(supp_top, funcs);
    }
  }
  else
  {
    static_assert( kMaxNumOutputs > 2, "Not implemented for more than 2 outputs" );
  }

  return std::nullopt;
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
auto SpinDecSolo<Ntk, Manager, Simulator, Params>::ChooseVars(std::vector<typename List::Signal>const& supp0, Reminders const& funcs, size_t num_vars) -> math::static_vector<typename List::Signal, Params::kMaxInpsSolo>
{
  std::vector<Function const *> funcs_bits;
  std::vector<Function const *> funcs_care;
  for ( auto const& func : funcs )
  {
    funcs_bits.push_back( &func._bits );
    funcs_care.push_back( &func._care );
  }
  spfds_.init( funcs_bits, funcs_care );
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
auto SpinDecSolo<Ntk, Manager, Simulator, Params>::ExtractReminders(std::vector<typename List::Signal> const& support) -> std::optional<RemindersType>
{
  spins_.clear();
  auto const& spfd = spfds_.Get(0);
  for ( uint32_t m = 0; m < spfd.num_masks_; ++m )
  {
    if (!spfd.is_killed(m))
    {
      Function const func = learner_.Run( spfd.func_[1], spfd.masks_[m], support );
      spins_.AddSpin( func, spfd.masks_[m], m );
    }
  }

  return GetFunctions();
}


  template <class Ntk, typename Manager, typename Simulator, typename Params>
  auto SpinDecSolo<Ntk, Manager, Simulator, Params>::GetFunctions() -> RemindersType
  {
    RemindersType res;
    if constexpr ( kMaxNumOutputs == 1u )
    {
      res[0] = Reminders();
      res[1] = Reminders();
    }
    else
    {
      res = Reminders();
    }

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

}  // namespace lut
}  // namespace synthesis
}  // namespace rinox

#endif