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

#include "../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../boolean/smoother.hpp"
#include "../../boolean/spfd.hpp"
#include "../../boolean/support_minimizer.hpp"
#include "../../math/ising.hpp"
#include "../../windowing/window_simulator.hpp"
#include "utils.hpp"

namespace rinox {
namespace synthesis {
namespace lut {

template <typename Params = DefaultLutSynthParams>
class SpinDecomposition {
  using CSFunction = kitty::static_truth_table<Params::kMaxCutsSize>;
  using ISFunction = kitty::ternary_truth_table<CSFunction>;
  using Reminder = math::static_vector<ISFunction, Params::kMaxNumOutputs>;
  using List =
      synthesis::lut::LutList<Params::kMaxInpsSolo, Params::kMaxNumOutputs,
                              Params::kMaxCutsSize>;
  using ListSignal = typename List::Signal;
  using ListSignals = math::static_vector<ListSignal, Params::kMaxNumOutputs>;
  using SuppMinimizer =
      rinox::boolean::support_minimizer<Params::kMaxCutsSize,
                                        Params::kExactSuppMin>;
  using Spfds = boolean::spfds_manager<CSFunction, 1u << Params::kMaxInpsSolo,
                                       Params::kMaxNumOutputs>;
  using Smoother = boolean::Smoother<CSFunction, typename List::Signal,
                                     Params::kMaxCutsSize>;

  static constexpr uint32_t kMaxNumSpins = 1u << Params::kMaxInpsSolo;

  using Spins = math::SpinSystem<CSFunction, kMaxNumSpins, 2u>;

 public:
  SpinDecomposition();
  bool Run(ISFunction const&);
  bool Run(Reminder const&);
  std::optional<ListSignals> RunRecursion(Reminder const&,
                                          std::vector<ListSignal>&);
  template <typename Ntk>
  std::vector<typename Ntk::signal> Insert(
      Ntk& ntk, std::vector<typename Ntk::signal> const& cut) {
    return list_.template Insert<Ntk>(ntk, cut);
  }

 private:
  List list_;
  SuppMinimizer supp_minimizer_;
  Spfds spfds_;
  Spins spins_;
  Smoother smoother_;
};

template <typename Params>
SpinDecomposition<Params>::SpinDecomposition() {
  auto const& proj_fns = ProjCache<Params::kMaxCutsSize>();
  for (int i = 0; i < Params::kMaxCutsSize; ++i) {
    list_.CreateInput(&proj_fns[i]);
  }
}

template <typename Params>
bool SpinDecomposition<Params>::Run(ISFunction const& func) {
  Reminder reminder;
  reminder.push_back(func);
  return Run(reminder);
}

template <typename Params>
bool SpinDecomposition<Params>::Run(Reminder const& reminder) {
  list_.RemoveGates();
  auto support = list_.GetInputSignals();
  auto res = RunRecursion(reminder, support);
  if (!res) return false;
  for (int i = 0; i < res->size(); ++i) list_.CreateOutput((*res)[i]);
  return true;
}

template <typename Params>
auto SpinDecomposition<Params>::RunRecursion(Reminder const& reminder,
                                             std::vector<ListSignal>& support)
    -> std::optional<ListSignals> {
  auto reminder_copy = reminder;
  std::vector<uint32_t> vars;
  for (auto const s : support) {
    vars.push_back(list_.GetIndex(s) - 1);
  }
  supp_minimizer_.template Run<Params::kMaxNumOutputs>(reminder_copy, vars,
                                                       support);
  if (support.size() <= Params::kMaxInpsSolo) {
    math::static_vector<ListSignal, Params::kMaxInpsSolo> fanin;
    for (auto const s : support) fanin.push_back(s);
    return list_.CreateInternal(fanin, reminder);
  }
  return std::nullopt;
}

}  // namespace lut
}  // namespace synthesis
}  // namespace rinox

#endif