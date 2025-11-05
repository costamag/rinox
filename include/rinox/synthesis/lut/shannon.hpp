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

#ifndef RINOX_SYNTHESIS_LUT_SHANNON_SOLO_HPP_
#define RINOX_SYNTHESIS_LUT_SHANNON_SOLO_HPP_

#include <fmt/format.h>

#include <chrono>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <optional>

#include "utils.hpp"
#include "../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../windowing/window_manager.hpp"
#include "../../windowing/window_simulator.hpp"

namespace rinox {
namespace synthesis {
namespace lut {

template <class Ntk, typename Manager, typename Simulator,
          typename Params = DefaultLutSynthParams>
class ShannonDecSolo {
  static constexpr uint32_t kMaxNumOutputs = Ntk::max_num_outputs;
  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using Func = kitty::static_truth_table<Params::kNumVarsCube>;
  using IFunc = kitty::ternary_truth_table<Func>;
  using Arrivals = analyzers::trackers::arrival_times_tracker<Ntk>;
  using List =
      LutList<Params::kMaxInpsSolo, kMaxNumOutputs, Params::kNumVarsCube>;

 public:
  ShannonDecSolo(Ntk& ntk, Manager& manager, Simulator& simulator, LutSynthStats& st)
      : ntk_(ntk), manager_(manager), simulator_(simulator), st_(st) {}

  template<uint32_t kNumInps>
  bool Run(IFunc const&, math::static_vector<Signal, kNumInps> const&);
  template<uint32_t kNumInps>
  bool Run(IFunc const&, math::static_vector<Signal, kNumInps> const&, Arrivals const&);
  std::optional<typename List::Signal> Decompose( IFunc, std::vector<uint8_t>);
  
  template<uint32_t kNumInps>
  Signal Insert(math::static_vector<Signal, kNumInps> const&);
  template<uint32_t kNumInps>
  void Init(math::static_vector<Signal, kNumInps> const&);
  template<uint32_t kNumInps>
  void Init(math::static_vector<Signal, kNumInps> const&, Arrivals const&);


  std::tuple<uint8_t, uint8_t> ChooseVar(std::vector<uint8_t> const&);

  size_t NumGates() const{ return list_.NumGates(); }
  double GetTime() const;
 private:
  Ntk& ntk_;
  Manager& manager_;
  Simulator& simulator_;
  LutSynthStats& st_;
  List list_;
};

// Implementation of Run without timing information
template <class Ntk, typename Manager, typename Simulator, typename Params>
template<uint32_t kNumInps>
bool ShannonDecSolo<Ntk, Manager, Simulator, Params>::Run(
    IFunc const& func, math::static_vector<Signal, kNumInps> const& cut) {
  bool success = mockturtle::call_with_stopwatch(st_.time, [&]() {
    Init<kNumInps>(cut);
    std::vector<uint8_t> vars;
    auto const num_vars = std::min(Params::kNumVarsCube, static_cast<uint32_t>(cut.size()));
    for (uint8_t v = 0u; v < num_vars; ++v) vars.push_back(v);
    auto s = Decompose(func, vars);
    if (s) {
      list_.CreateOutput(*s);
      return true;
    } else
      return false;
  });
  return success;
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
template<uint32_t kNumInps>
bool ShannonDecSolo<Ntk, Manager, Simulator, Params>::Run(
    IFunc const& func, math::static_vector<Signal, kNumInps> const& cut,
    Arrivals const& arrivals) {
  bool success = mockturtle::call_with_stopwatch(st_.time, [&]() {
    Init<kNumInps>(cut, arrivals);
    std::vector<uint8_t> vars;
    auto const num_vars = std::min(Params::kNumVarsCube, static_cast<uint32_t>(cut.size()));
    for (uint8_t v = 0u; v < num_vars; ++v) vars.push_back(v);
    auto s = Decompose(func, vars);
//kitty::print_binary(func); std::cout << std::endl;
    if (s) {
      list_.CreateOutput(*s);
      return true;
    } else
      return false;
  });
  return success;
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
template<uint32_t kNumInps>
void ShannonDecSolo<Ntk, Manager, Simulator, Params>::Init(
    math::static_vector<Signal, kNumInps> const& cut) {
  list_.Reset(manager_.num_mffc());
  // assign the input signatures to the LUT list
  for (Signal const f : cut) list_.CreateInput(&simulator_.Get(f));
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
template<uint32_t kNumInps>
void ShannonDecSolo<Ntk, Manager, Simulator, Params>::Init(
    math::static_vector<Signal, kNumInps> const& cut, Arrivals const& arrivals) {
  list_.Reset( manager_.num_mffc() );
  // assign the input signatures to the LUT list
  for (Signal const f : cut)
    list_.CreateInput(&simulator_.Get(f), arrivals.get_time(f));
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
std::optional<typename ShannonDecSolo<Ntk, Manager, Simulator, Params>::List::Signal>
ShannonDecSolo<Ntk, Manager, Simulator, Params>::Decompose(
    IFunc func,
    std::vector<uint8_t> vars) {
  if ( list_.NumGates() > manager_.num_mffc() ) return std::nullopt;
  if (vars.size() <= Params::kMaxInpsSolo) {
    math::static_vector<IFunc, kMaxNumOutputs> funcs;
    funcs.push_back( func );

    math::static_vector<typename List::Signal, Params::kMaxInpsSolo> fanins;
    for (uint8_t v : vars)
    {
      fanins.push_back(list_.MakeSignal(static_cast<typename List::Index>(v), 0));
    }

    return list_.CreateInternal(fanins, funcs);
  }

  auto [i, v] = ChooseVar(vars);
  vars.erase(vars.begin() + i);
  
  typename List::Signal x = list_.MakeSignal(static_cast<typename List::Index>(v), 0);
  Func const& sim_x = list_.GetSimBits(x);
  IFunc func1 = {func._bits,func._care & sim_x};
  IFunc func0 = {func._bits,func._care & ~sim_x};
  auto f1 = Decompose(func1, vars);
  auto f0 = Decompose(func0, vars);
  if (!f0 || !f1) return std::nullopt;
  Func const& sim_1 = list_.GetSimBits(*f1);
  Func const& sim_0 = list_.GetSimBits(*f0);
  Func bits = ((sim_x & sim_1) | ((~sim_x) & sim_0));
  Func care = func._care;

  //assert( kitty::equal( sim & func._care, func._bits & func._care ) );
  math::static_vector<typename List::Signal, Params::kMaxInpsSolo> fanins;
  fanins.push_back(x) ;
  fanins.push_back(*f1);
  fanins.push_back(*f0);
  math::static_vector<IFunc, kMaxNumOutputs> funcs;
  funcs.push_back( IFunc{bits, care} );

  return list_.CreateInternal(fanins, funcs);
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
std::tuple<uint8_t, uint8_t>
ShannonDecSolo<Ntk, Manager, Simulator, Params>::ChooseVar(
    std::vector<uint8_t> const& vars) {
  uint8_t v = 0;
  double latest = std::numeric_limits<double>::min();
  for (uint8_t i = 0; i < vars.size(); ++i) {
    double const time =
        list_.NodeTime(static_cast<typename List::Index>(vars[i]));
    if (time > latest) {
      latest = time;
      v = i;
    }
  }
  return std::tuple{v, vars[v]};
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
template<uint32_t kNumInputs>
typename Ntk::signal ShannonDecSolo<Ntk, Manager, Simulator, Params>::Insert(
    math::static_vector<Signal, kNumInputs> const& cut) {
//std::cout << "inserting: ";
//for ( auto f : cut )
//  std::cout << f.index << "." << f.output << " ";
//std::cout << std::endl;

  std::unordered_map<typename List::Signal, typename Ntk::signal> map;
  for ( uint32_t i = 0; i < cut.size(); ++i )
    map[list_.MakeSignal(i,0)] = cut[i];

  list_.ForeachGate( [&]( auto const& internal_node, auto index ){
    std::visit([&](auto const& n) {
      using T = std::decay_t<decltype(n)>;
      if constexpr (std::is_same_v<T, typename List::Internal>) {
//std::cout << index << ": ";
        auto const& fanins = n.fanins_;
        auto const& sims = n.sims_;
        std::vector<typename List::Func const *> ptrs;
        std::vector<Signal> children;
        for ( auto f : fanins )
        {
          children.push_back( map[f] );
          ptrs.push_back( &list_.GetSimBits(f) );
        }
        auto tts = Interpolate<Func, kMaxNumOutputs>( ptrs, sims );
        auto fnew = ntk_.create_node( children, tts );
//for ( auto f : children )
//  std::cout << f.index << "." << f.output << " ";
//std::cout << std::endl;
        for ( uint32_t i=0; i<tts.size(); ++i )
        {
//kitty::print_binary(tts[i]);std::cout << " " << index << "." << i << std::endl;
          auto s = list_.MakeSignal(index,i);
          typename Ntk::signal f{ fnew.index, i };
          map[s] = f;
          auto const& sim = ntk_.compute( f, ptrs );
          list_.SetSimBits( s, sim );
          list_.SetSimCare( s, sim|~sim);
        }
      }
    }, internal_node);
  } );
  auto outputs = list_.GetOutputs();
  assert( outputs.size() == 1u );
  return map[outputs[0]];
}

template <class Ntk, typename Manager, typename Simulator, typename Params>
double ShannonDecSolo<Ntk, Manager, Simulator, Params>::GetTime() const{
  return list_.GetTime();
}

}  // namespace lut
}  // namespace synthesis
}  // namespace rinox

#endif