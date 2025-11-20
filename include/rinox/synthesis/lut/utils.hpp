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

template <typename TT, uint32_t NumVars>
std::vector<kitty::dynamic_truth_table> Interpolate(
    std::vector<TT const*> const& tts,
    math::static_vector<TT, NumVars> const& target) {
  std::vector<kitty::dynamic_truth_table> res;
  for (int i = 0; i < target.size(); ++i) {
    const uint64_t num_vars = tts.size();
    kitty::dynamic_truth_table func(num_vars);
    TT tmp;
    for (uint64_t m = 0; m < (1u << num_vars); ++m) {
      tmp |= ~tmp;
      for (uint32_t v = 0; v < num_vars; ++v) {
        tmp &= (((m >> v) & 0x1) == 0x1) ? *tts[v] : ~(*tts[v]);
      }
      if (kitty::count_ones(target[i] & tmp) > 0) kitty::set_bit(func, m);
    }
    res.push_back(func);
  }
  return res;
}

template <uint32_t N>
inline const std::array<kitty::static_truth_table<N>, N>& ProjCache() {
  using TTF = kitty::static_truth_table<N>;
  static const std::array<TTF, N> cache = [] {
    std::array<TTF, N> a{};
    for (uint32_t i = 0; i < N; ++i) kitty::create_nth_var(a[i], i);
    return a;
  }();
  return cache;
}

template <uint32_t NumVarsFunction, uint32_t NumVarsSignature,
          uint32_t kMaxNumOutputs>
math::static_vector<
    kitty::ternary_truth_table<kitty::static_truth_table<NumVarsFunction>>,
    kMaxNumOutputs>
Interpolate(
    std::vector<const kitty::static_truth_table<NumVarsSignature>*> const& tts,
    math::static_vector<kitty::static_truth_table<NumVarsSignature>,
                        kMaxNumOutputs> const& onset,
    math::static_vector<kitty::static_truth_table<NumVarsSignature>,
                        kMaxNumOutputs> const& careset) {
  using TTS = kitty::static_truth_table<NumVarsSignature>;
  using TTF = kitty::static_truth_table<NumVarsFunction>;

  auto const& proj_funcs = ProjCache<NumVarsFunction>();  // cached once per N
  math::static_vector<
      kitty::ternary_truth_table<kitty::static_truth_table<NumVarsFunction>>,
      kMaxNumOutputs>
      res;

  for (uint32_t i = 0; i < onset.size(); ++i) {
    TTF func;
    TTF care;
    TTF mask_function;
    TTS mask_signature;
    const uint64_t num_vars = tts.size();
    for (uint64_t m = 0; m < (1u << num_vars); ++m) {
      mask_function |= ~mask_function;
      mask_signature |= ~mask_signature;
      for (uint32_t v = 0; v < num_vars; ++v) {
        mask_signature &= (((m >> v) & 0x1) == 0x1) ? *tts[v] : ~(*tts[v]);
        mask_function &=
            (((m >> v) & 0x1) == 0x1) ? proj_funcs[v] : ~proj_funcs[v];
      }
      if (kitty::count_ones(careset[i] & mask_signature) > 0)
        care |= mask_function;
      if (kitty::count_ones(onset[i] & mask_signature) > 0)
        func |= mask_function;
    }
    kitty::ternary_truth_table<kitty::static_truth_table<NumVarsFunction>>
        result(func, care);
    res.push_back(result);
  }
  return res;
}

template <typename TT, uint32_t NumVars>
std::vector<kitty::dynamic_truth_table> Interpolate(
    std::vector<TT const*> const& tts,
    math::static_vector<kitty::ternary_truth_table<TT>, NumVars> const&
        target) {
  std::vector<kitty::dynamic_truth_table> res;
  for (int i = 0; i < target.size(); ++i) {
    const uint64_t num_vars = tts.size();
    kitty::dynamic_truth_table func(num_vars);
    TT tmp;
    for (uint64_t m = 0; m < (1u << num_vars); ++m) {
      tmp |= ~tmp;
      for (uint32_t v = 0; v < num_vars; ++v) {
        tmp &= (((m >> v) & 0x1) == 0x1) ? *tts[v] : ~(*tts[v]);
      }
      auto const mask = target[i]._care & tmp;
      auto const num_0s = kitty::count_ones(~target[i]._bits & mask);
      auto const num_1s = kitty::count_ones(target[i]._bits & mask);
      if (num_1s == 0)
        kitty::clear_bit(func, m);
      else if (num_0s == 0)
        kitty::set_bit(func, m);
      else
        std::cout << num_1s << " " << num_0s << " both non zero" << std::endl;
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
  static constexpr uint32_t kMaxCutsSize = 8u;
  static constexpr uint32_t kMaxInpsSolo = 6u;
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr uint32_t kMaxNumGates = 32u;
  static constexpr uint32_t kMaxNumOutputs = 2u;
  static constexpr bool kExactSuppMin = true;
};

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
class LutList {
 public:
  using Func = kitty::static_truth_table<kNumVarsCube>;
  using IFunc = kitty::ternary_truth_table<Func>;
  using Signal = uint64_t;
  using Index = uint64_t;
  static constexpr uint32_t kOutputBits = math::log2_ceil(kMaxNumOuts);
  static constexpr uint32_t kOutputMask = (1u << kOutputBits) - 1;

  struct Const {
    bool value_;
    double time_{0.0};
  };

  struct Input {
    double time_;
    Func const* sim_;
  };

  struct Internal {
    double time_;
    math::static_vector<Signal, kFaninSize> fanins_;
    math::static_vector<IFunc, kMaxNumOuts> sims_;
  };

  using LutNode = std::variant<Input, Internal, Const>;

  LutList();
  void Reset(size_t);
  void RemoveGates();
  Index CreateInput(Func const*, double time = 0.0);
  math::static_vector<Signal, kMaxNumOuts> CreateInternal(
      math::static_vector<Signal, kFaninSize> const&,
      math::static_vector<IFunc, kMaxNumOuts> const&);
  void CreateOutput(Signal);
  std::vector<Signal> GetInputSignals() const { return inputs_; };
  std::vector<Signal> GetInputIndices() const {
    std::vector<Index> inputs;
    for (auto const s : inputs_) inputs.push_back(GetIndex(s) - 1);
    return inputs;
  };
  Func const& GetSimBits(Signal) const;
  void SetSimBits(Signal, Func const&);
  void SetSimCare(Signal, Func const&);
  std::vector<Signal> GetOutputs() const { return outputs_; };

  template <typename Ntk>
  std::vector<typename Ntk::signal> Insert(
      Ntk& ntk, std::vector<typename Ntk::signal> const&);

  template <typename Fn>
  void ForeachGate(Fn&&);

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
  bool IsConst(Index i) const {
    return std::holds_alternative<Const>(nodes_[i]);
  }
  bool IsInternal(Index i) const {
    return std::holds_alternative<Internal>(nodes_[i]);
  }
  double NodeTime(Index i) const;
  double GetTime() const;

  size_t NumGates() const { return num_gates_; }

 private:
  uint32_t num_gates_ = 0;
  std::vector<LutNode> nodes_;
  std::vector<Signal> outputs_;
  std::vector<Signal> inputs_;
};

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::LutList() {
  Reset(0);
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
void LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::Reset(size_t size) {
  nodes_.clear();
  num_gates_ = 0;
  nodes_.reserve(size + 1);
  nodes_.push_back(Const{false});
  outputs_.clear();
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
void LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::RemoveGates() {
  num_gates_ = 0;
  outputs_.clear();
  for (int i = nodes_.size() - 1; i >= 0; --i) {
    if (std::is_same_v<std::decay_t<decltype(nodes_[i])>, Internal>)
      nodes_.erase(nodes_.begin() + i);
  }
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
typename LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::Index
LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::CreateInput(Func const* sim,
                                                            double time) {
  auto const index = static_cast<Index>(nodes_.size());
  nodes_.push_back(Input{time, sim});
  inputs_.push_back(MakeSignal(index, 0));
  return index;
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
auto LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::CreateInternal(
    math::static_vector<Signal, kFaninSize> const& fanins,
    math::static_vector<IFunc, kMaxNumOuts> const& sims)
    -> math::static_vector<Signal, kMaxNumOuts> {
  num_gates_++;
  double latest = 0;
  for (auto f : fanins) latest = std::max(NodeTime(GetIndex(f)), latest);

  auto const index = static_cast<Index>(nodes_.size());
  nodes_.push_back(Internal{latest + 1, fanins, sims});
  math::static_vector<Signal, kMaxNumOuts> res;
  for (uint32_t out = 0; out < sims.size(); ++out)
    res.push_back(MakeSignal(index, out));
  return res;
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
  for (auto const s : outputs_) {
    const auto idx = GetIndex(s);
    const auto out = GetOutput(s);
    double time = std::visit(
        [&](auto const& n) -> double {
          if constexpr (std::is_same_v<std::decay_t<decltype(n)>, Input>) {
            return n.time_;  // PI: external sim
          } else {
            return n.time_;  // Internal: owned sim
          }
        },
        nodes_[idx]);
    latest = std::min(latest, time);
  }
  return latest;
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
typename LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::Func const&
LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::GetSimBits(Signal s) const {
  const auto idx = GetIndex(s);
  const auto out = GetOutput(s);
  return std::visit(
      [&](auto const& n) -> Func const& {
        if constexpr (std::is_same_v<std::decay_t<decltype(n)>, Input>) {
          return *n.sim_;  // PI: external sim
        } else if constexpr (std::is_same_v<std::decay_t<decltype(n)>,
                                            Internal>) {
          return n.sims_[out]._bits;  // Internal: owned sim
        } else {
          throw std::runtime_error("[e] Trying to access const simulation");
        }
      },
      nodes_[idx]);
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
void LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::SetSimBits(
    Signal s, Func const& func) {
  const auto idx = GetIndex(s);
  const auto out = GetOutput(s);
  std::visit(
      [&](auto& n) -> void {
        if constexpr (std::is_same_v<std::decay_t<decltype(n)>, Input>) {
          throw std::runtime_error("[e] Trying to assign input simulation");
        } else {
          n.sims_[out]._bits = func;  // Internal: owned sim
        }
      },
      nodes_[idx]);
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
void LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::SetSimCare(
    Signal s, Func const& func) {
  const auto idx = GetIndex(s);
  const auto out = GetOutput(s);
  std::visit(
      [&](auto& n) -> void {
        if constexpr (std::is_same_v<std::decay_t<decltype(n)>, Input>) {
          throw std::runtime_error("[e] Trying to assign input simulation");
        } else {
          n.sims_[out]._care = func;  // Internal: owned sim
        }
      },
      nodes_[idx]);
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
template <typename Fn>
void LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::ForeachGate(Fn&& fn) {
  for (uint32_t index = 0u; index < nodes_.size(); ++index) {
    auto& n = nodes_[index];
    fn(n, index);
  }
}

template <uint32_t kFaninSize, uint32_t kMaxNumOuts, uint32_t kNumVarsCube>
template <typename Ntk>
std::vector<typename Ntk::signal>
LutList<kFaninSize, kMaxNumOuts, kNumVarsCube>::Insert(
    Ntk& ntk, std::vector<typename Ntk::signal> const& cut) {
  std::vector<typename Ntk::signal> res;
  std::unordered_map<Signal, typename Ntk::signal> index_to_signal;
  for (int i = 0; i < cut.size(); ++i) {
    index_to_signal[MakeSignal(i + 1, 0)] = cut[i];
  }

  ForeachGate([&](LutNode& v, std::size_t i) {
    if (auto* p = std::get_if<Internal>(&v)) {
      auto& n = *p;

      std::vector<typename Ntk::signal> fanins;
      fanins.reserve(n.fanins_.size());

      std::vector<Func const*> sims_ptrs;
      sims_ptrs.reserve(n.fanins_.size());

      for (auto f : n.fanins_) {
        fanins.push_back(index_to_signal.at(f));
        sims_ptrs.push_back(&GetSimBits(f));
      }

      const auto tts = Interpolate<Func, kMaxNumOuts>(sims_ptrs, n.sims_);
      const auto new_signals = ntk.create_node(fanins, tts);

      for (std::size_t out = 0; out < tts.size(); ++out) {
        index_to_signal[MakeSignal(i, out)] =
            typename Ntk::signal{new_signals.index, static_cast<unsigned>(out)};
      }

      auto const sims = ntk.compute(new_signals.index, sims_ptrs);
      for (int i = 0; i < sims.size(); ++i) {
        n.sims_[i]._bits = sims[i];
        n.sims_[i]._care |= ~n.sims_[i]._care;
      }
    }
  });

  for (auto const s : outputs_) {
    assert(index_to_signal.find(s) != index_to_signal.end() &&
           "output signal not found");
    res.push_back(index_to_signal.at(s));
  }
  return res;
}

template <typename StaticParams = DefaultLutSynthParams>
struct SubProblems {
  using Signature = kitty::static_truth_table<StaticParams::kMaxCutsSize>;
  using Function = kitty::ternary_truth_table<Signature>;

  void Insert(std::vector<uint32_t> const& cut, Function const& func) {
    cuts.push_back(cut);
    funcs.push_back({func});
  }

  void clear() {
    cuts.clear();
    funcs.clear();
  }

  math::static_vector<std::vector<uint32_t>, StaticParams::kMaxNumGates> cuts;
  math::static_vector<
      math::static_vector<Function, StaticParams::kMaxNumOutputs>,
      StaticParams::kMaxNumGates>
      funcs;
};

}  // namespace lut
}  // namespace synthesis
}  // namespace rinox

#endif