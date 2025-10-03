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
  \file lut_list.hpp
  \brief Index list for LUT networks

  \author Andrea Costamagna
*/

#ifndef RINOX_SYNTHESIS_LUT_LIST_HPP_
#define RINOX_SYNTHESIS_LUT_LIST_HPP_

#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <optional>
#include <variant>

#include "../../math/math.hpp"

namespace rinox {
namespace synthesis {
namespace lut {

template<typename TT, uint32_t NumVars>
std::vector<kitty::dynamic_truth_table> Interpolate( std::vector<TT const*> const& tts, math::static_vector<TT, NumVars> const& target )
{
  std::vector<kitty::dynamic_truth_table> res;
  for ( int i = 0; i < target.size(); ++i )
  {
    const uint64_t num_vars = tts.size();
    kitty::dynamic_truth_table func(num_vars);
    TT tmp;
    for ( uint64_t m = 0; m < (1u<<num_vars); ++m )
    {
      tmp |= ~tmp;
      for ( uint32_t v=0; v<num_vars; ++v )
      {
        tmp &= ( ( ( m>>v ) & 0x1 ) == 0x1 ) ? *tts[v] : ~( *tts[v] );
      }
      if ( kitty::count_ones( target[i]&tmp ) > 0 )
        kitty::set_bit( func, m );
    }
    res.push_back(func);
  }
  return res;
}

template<typename TT, uint32_t NumVars>
std::vector<kitty::dynamic_truth_table> Interpolate( std::vector<TT const*> const& tts, math::static_vector<kitty::ternary_truth_table<TT>, NumVars> const& target )
{
  std::vector<kitty::dynamic_truth_table> res;
  for ( int i = 0; i < target.size(); ++i )
  {
    const uint64_t num_vars = tts.size();
    kitty::dynamic_truth_table func(num_vars);
    TT tmp;
    for ( uint64_t m = 0; m < (1u<<num_vars); ++m )
    {
      tmp |= ~tmp;
      for ( uint32_t v=0; v<num_vars; ++v )
      {
        tmp &= ( ( ( m>>v ) & 0x1 ) == 0x1 ) ? *tts[v] : ~( *tts[v] );
      }
      if ( kitty::count_ones( target[i]._bits&target[i]._care&tmp ) > 0 )
        kitty::set_bit( func, m );
    }
    res.push_back(func);
  }
  return res;
}

struct LutSynthStats {
  mockturtle::stopwatch<>::duration time{0};

  double GetTime() const { return mockturtle::to_seconds(time); }

  void Report() const {
    std::cout << fmt::format("============ SHANNON DEC =============\n");
    std::cout << fmt::format("    time                          = {:5.2f}\n",
                             GetTime());
  }
};

struct DefaultLutSynthParams {
  static constexpr uint32_t kNumVarsCube = 12u;
  static constexpr uint32_t kMaxInpsSolo = 6u;
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr uint32_t kMaxNumGates = 32u;
};

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
class LutList {
 public:
  using Func = kitty::static_truth_table<kNumVarsCube>;
  using IFunc = kitty::ternary_truth_table<Func>;
  using Signal = uint32_t;
  using Index = uint32_t;
  static constexpr uint32_t kOutputBits = math::log2_ceil(kMaxNumOuts);
  static constexpr uint32_t kOutputMask = (1u << kOutputBits) - 1;

  struct Input {
    double time_;
    Func const* sim_;
  };

  struct Internal {
    double time_;
    math::static_vector<Signal, kFaninSize> fanins_;
    math::static_vector<IFunc, kMaxNumOuts> sims_;
  };

  using LutNode = std::variant<Input, Internal>;

  void Reset(size_t);
  Index CreateInput(Func const*, double time = 0.0);
  Signal CreateInternal(math::static_vector<Signal, kFaninSize> const&,
                       math::static_vector<IFunc, kMaxNumOuts> const&);
  void CreateOutput(Signal);

  Func const& GetSimBits(Signal) const;
  void SetSimBits(Signal, Func const&);
  void SetSimCare(Signal, Func const&);
  std::vector<Signal> GetOutputs() const { return outputs_; };

  template<typename Fn>
  void ForeachGate( Fn&& ) const;

  static constexpr Index GetIndex(Signal f) {
    return static_cast<Index>(f >> kOutputBits);
  }
  static constexpr Index GetOutput(Signal f) {
    return static_cast<Index>(f & kOutputMask);
  }
  static constexpr Signal MakeSignal(Index idx, uint32_t out) {
    return (static_cast<Signal>(idx) << kOutputBits) | (out & kOutputMask);
  }

  bool IsInput(Index i) const {
    return std::holds_alternative<Input>(nodes_[i]);
  }
  double NodeTime(Index i) const;
  double GetTime() const;

  size_t NumGates() const { return num_gates_; }

 private:
  uint32_t num_gates_ = 0;
  std::vector<LutNode> nodes_;
  std::vector<Signal> outputs_;
};

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
void LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::Reset(size_t size) {
  num_gates_ = 0;
  nodes_.reserve(size);
  nodes_.clear();
  outputs_.clear();
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
typename LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::Index LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::CreateInput(
    Func const* sim, double time) {
  auto const index = static_cast<Index>(nodes_.size());
  nodes_.push_back(Input{time, sim});
  return index;
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
auto LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::CreateInternal(
    math::static_vector<Signal, kFaninSize> const& fanins,
    math::static_vector<IFunc, kMaxNumOuts> const& sims) -> Signal {
  
  num_gates_++;
  double latest = 0;
  for (auto f : fanins) latest = std::max(NodeTime(GetIndex(f)), latest);

  auto const index = static_cast<Index>(nodes_.size());
  nodes_.push_back(Internal{latest + 1, fanins, sims});
  return MakeSignal(index,0);
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
void LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::CreateOutput(Signal f) {
  outputs_.emplace_back(f);
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
double LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::NodeTime(Index i) const {
  return std::visit([](auto const& n) { return n.time_; }, nodes_[i]);
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
double LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::GetTime() const {
  double latest = std::numeric_limits<double>::min();
  for ( auto const s : outputs_ )
  {
    const auto idx = GetIndex(s);
    const auto out = GetOutput(s);
    double time = std::visit([&](auto const& n) -> double {
      if constexpr (std::is_same_v<std::decay_t<decltype(n)>, Input>) {
        return n.time_; // PI: external sim
      } else {
        return n.time_; // Internal: owned sim
      }
    }, nodes_[idx]);
    latest = std::min( latest, time );
  }
  return latest;
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
typename LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::Func const& LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::GetSimBits(Signal s) const {
  const auto idx = GetIndex(s);
  const auto out = GetOutput(s);
  return std::visit([&](auto const& n) -> Func const& {
    if constexpr (std::is_same_v<std::decay_t<decltype(n)>, Input>) {
      return *n.sim_; // PI: external sim
    } else {
      return n.sims_[out]._bits; // Internal: owned sim
    }
  }, nodes_[idx]);
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
void LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::SetSimBits(Signal s, Func const& func )
{
  const auto idx = GetIndex(s);
  const auto out = GetOutput(s);
  std::visit([&](auto & n) -> void {
    if constexpr (std::is_same_v<std::decay_t<decltype(n)>, Input>) {
      throw std::runtime_error("[e] Trying to assign input simulation");
    } else {
      n.sims_[out]._bits = func; // Internal: owned sim
    }
  }, nodes_[idx]);
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
void LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::SetSimCare(Signal s, Func const& func )
{
  const auto idx = GetIndex(s);
  const auto out = GetOutput(s);
  std::visit([&](auto & n) -> void {
    if constexpr (std::is_same_v<std::decay_t<decltype(n)>, Input>) {
      throw std::runtime_error("[e] Trying to assign input simulation");
    } else {
      n.sims_[out]._care = func; // Internal: owned sim
    }
  }, nodes_[idx]);
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
template<typename Fn>
void LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::ForeachGate(Fn && fn) const {
  for ( uint32_t index = 0u; index < nodes_.size(); ++index )
  {
    auto const& n = nodes_[index];
    if ( !IsInput( index ) )
      fn( n, index );
  }
}

}  // namespace lut
}  // namespace synthesis
}  // namespace rinox

#endif