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
  math::static_vector<ListSignal, Params::kMaxInpsSolo> ChooseVars(
      std::vector<ListSignal> const&, Reminder const&, uint32_t);
  template <typename Ntk>
  std::vector<typename Ntk::signal> Insert(
      Ntk& ntk, std::vector<typename Ntk::signal> const& cut) {
    return list_.template Insert<Ntk>(ntk, cut);
  }

  bool SimplerThan(std::vector<ListSignal> const& supp_a, Reminder const& a,
                   std::vector<ListSignal> const& supp_b, Reminder const& b) {
    if (supp_a.size() < supp_b.size()) return true;
    uint32_t cost_a = 0;
    uint32_t cost_b = 0;
    for (auto const& f : a) cost_a += kitty::count_ones(f._care);
    for (auto const& f : b) cost_b += kitty::count_ones(f._care);
    return cost_a < cost_b;
  }

  Reminder ExtractSimplest();
  Reminder ExtractSimplest2();

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
  auto vars_rem = vars;
  auto support_rem = support;
  supp_minimizer_.template Run<Params::kMaxNumOutputs>(reminder_copy, vars,
                                                       support);
  if (support.size() <= Params::kMaxInpsSolo) {
    math::static_vector<ListSignal, Params::kMaxInpsSolo> fanin;
    for (auto const s : support) fanin.push_back(s);
    return list_.CreateInternal(fanin, reminder);
  }
  if constexpr (Params::kMaxNumOutputs == 1) {
    auto supp_top =
        ChooseVars(support, reminder_copy, Params::kMaxInpsSolo - 1);
    Reminder smooth_func = ExtractSimplest();

    supp_minimizer_.template Run<Params::kMaxNumOutputs>(smooth_func, vars_rem,
                                                         support_rem);

    if (SimplerThan(support_rem, smooth_func, support, reminder)) {
      std::vector<CSFunction const*> sims_ptrs;
      for (int i = 0; i < smooth_func.size(); ++i)
        sims_ptrs.push_back(&smooth_func[i]._bits);
      auto f = RunRecursion(smooth_func, support_rem);
      if (f) {
        for (auto s : *f) {
          supp_top.push_back(s);
        }
        return list_.CreateInternal(supp_top, reminder);
      }
    } else {
      std::cout << "not simpler than" << std::endl;
    }
  }

  return std::nullopt;
}

template <typename Params>
auto SpinDecomposition<Params>::ChooseVars(
    std::vector<ListSignal> const& support, Reminder const& reminder,
    uint32_t num_vars)
    -> math::static_vector<ListSignal, Params::kMaxInpsSolo> {
  math::static_vector<ListSignal, Params::kMaxInpsSolo> res;
  std::vector<CSFunction const*> bits_vec;
  std::vector<CSFunction const*> care_vec;
  for (int i = 0; i < reminder.size(); ++i) {
    bits_vec.push_back(&reminder[i]._bits);
    care_vec.push_back(&reminder[i]._care);
  }
  spfds_.init(bits_vec, care_vec);
  while (!spfds_.is_covered() && res.size() < num_vars) {
    std::optional<ListSignal> best;
    auto best_cost = spfds_.get_num_edges();
    for (auto const& s : support) {
      auto const& sim = list_.GetSimBits(s);
      auto cost = spfds_.evaluate(sim);
      if (cost < best_cost) {
        best_cost = cost;
        best = s;
      }
    }
    if (best) {
      res.push_back(*best);
      spfds_.update(list_.GetSimBits(*best));
    }
  }
  return res;
}

template <typename Params>
auto SpinDecomposition<Params>::ExtractSimplest() -> Reminder {
  Reminder res;
  for (int i = 0; i < Params::kMaxNumOutputs; ++i) {
    spins_.clear();
    auto const& spfd = spfds_.Get(i);

    CSFunction care;
    for (uint32_t m = 0; m < spfd.get_num_masks(); ++m) {
      if (!spfd.is_killed(m)) care |= spfd.masks_[m];
    }

    for (uint32_t m = 0; m < spfd.get_num_masks(); ++m) {
      if (!spfd.is_killed(m)) {
        auto const info =
            smoother_.RunLabelPropagation(spfd.func_[1], spfd.masks_[m]);
        spins_.AddSpin(info & care, care, m);
      }
    }
    spins_.MinimizeEnergy();
    ISFunction func;
    spins_.ForeachSpin([&](auto const& spin) {
      if (spin.value > 0)
        func._bits |= spfd.masks_[spin.name] & spfd.func_[1];
      else
        func._bits |= spfd.masks_[spin.name] & spfd.func_[0];
    });
    func._care = care;
    res.push_back(func);
  }
  return res;
}

template <typename Params>
auto SpinDecomposition<Params>::ExtractSimplest2() -> Reminder {
  Reminder res;
  for (int i = 0; i < Params::kMaxNumOutputs; ++i) {
    spins_.clear();
    auto const& spfd = spfds_.Get(i);

    CSFunction care;
    for (uint32_t m = 0; m < spfd.get_num_masks(); ++m) {
      if (!spfd.is_killed(m)) care |= spfd.masks_[m];
    }

    for (uint32_t m = 0; m < spfd.get_num_masks(); ++m) {
      if (!spfd.is_killed(m)) {
        auto const info =
            smoother_.RunLabelPropagation(spfd.func_[1], spfd.masks_[m]);
        spins_.AddSpin(info & care, care, m);
      }
    }
    spins_.MinimizeEnergy();
    ISFunction func;
    spins_.ForeachSpin([&](auto const& spin) {
      if (spin.value > 0)
        func._bits |= spfd.masks_[spin.name] & spfd.func_[1];
      else
        func._bits |= spfd.masks_[spin.name] & spfd.func_[0];
    });
    func._care = care;
    res.push_back(func);
  }
  return res;
}

}  // namespace lut
}  // namespace synthesis
}  // namespace rinox

#endif