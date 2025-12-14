/* rinox: C++ logic network library
 * Copyright (C) 20X8-2023  EPFL
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

#ifndef RINOX_ALGORITHMS_CELL_REWRITE_DELAY_HPP_
#define RINOX_ALGORITHMS_CELL_REWRITE_DELAY_HPP_

#include <fmt/format.h>

#include <chrono>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <optional>

#include "../../../analyzers/evaluators/times_evaluator.hpp"
#include "../../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../../boolean/simd.hpp"
#include "../../../boolean/spfd.hpp"
#include "../../../dependency/struct_dependencies.hpp"
#include "../../../dependency/window/oracles.hpp"
#include "../../../network/visualize.hpp"
#include "../../../synthesis/lut/synth.hpp"
#include "../../../windowing/window_manager_tiered.hpp"
#include "../../../windowing/window_simulator_tiered.hpp"
#include "../../profilers/area_profiler.hpp"
#include "../../profilers/delay_profiler.hpp"
#include "../../profilers/power_profiler.hpp"
#include "../../profilers/profilers_utils.hpp"

#ifndef RINOX_NUM_VARS_SIGN
#define RINOX_NUM_VARS_SIGN 6
#endif

#ifndef RINOX_MAX_CUTS_SIZE
#define RINOX_MAX_CUTS_SIZE 6
#endif

#ifndef RINOX_MAX_CUBE_SPFD
#define RINOX_MAX_CUBE_SPFD 12
#endif

#ifndef RINOX_MAX_NUM_LEAVES
#define RINOX_MAX_NUM_LEAVES 12
#endif

namespace rinox {
namespace opto {
namespace cell {

struct RewriteDelayParamsDefault : windowing::WindowManagerParamsDefault {
  static constexpr uint32_t kMaxNumLeaves = 8;
  static constexpr uint32_t kInpsCapacity = 512u;
  static constexpr uint32_t kDivsCapacity = 100u;
  static constexpr uint32_t kMffcCapacity = 512u;
  static constexpr uint32_t kTfosCapacity = 512u;
  static constexpr uint32_t kOutsCapacity = 512u;
  static constexpr uint32_t kMaxCutsSize = 4u;
  static constexpr uint32_t kMaxCutReconv = 3u;
  static constexpr uint32_t kMaxNumVarsDb = 4u;
  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 0;
  int max_num_tiers_tfi = 10;
  int max_num_tiers_tfo = 0;
  bool area_constrained = false;
};

struct RewriteDelayStatsDefault : windowing::WindowManagerStatsDefault {
  mockturtle::stopwatch<>::duration time_rewrite{0};
  uint32_t num_rewrite{0};

  double GetTimeRewrite() const { return mockturtle::to_seconds(time_rewrite); }

  double GetTimeTotal() const { return GetTimeRewrite(); }

  void RecordRewrite() { num_rewrite++; }

  uint32_t GetNumRewrite() const { return num_rewrite; }

  void Report() const {
    std::cout << fmt::format("================ TIME  ====================\n");
    std::cout << fmt::format("    time rewrite                   = {:5.2f}\n",
                             GetTimeRewrite());
    std::cout << fmt::format("    time total                     = {:5.2f}\n",
                             GetTimeTotal());
    std::cout << fmt::format("================ QUALITY ==================\n");
    std::cout << fmt::format("    count rewrite                  = {:5d}\n",
                             GetNumRewrite());
  }
};

namespace detail {

template <class Ntk, typename Database, typename Params, typename Stats>
class RewriteDelayImpl {
  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using Signals = rinox::math::static_vector<Signal, Params::kMaxCutsSize>;
  using Signature = kitty::static_truth_table<Params::kMaxNumLeaves>;
  using WindowManagerT = windowing::WindowManager<Ntk, Params>;
  using WindowSimulatorT =
      windowing::WindowSimulator<Ntk, Params::kMaxNumLeaves>;
  static constexpr uint32_t kMaxNumOutputs = Ntk::kMaxNumOutputs;
  using TimesEvaluatorT = typename analyzers::evaluators::TimesEvaluator<Ntk>;
  using Chain = evaluation::chains::bound_chain<
      rinox::network::design_type_t::CELL_BASED>;
  using StructDependencies = dependency::struct_dependencies<Ntk, Params>;
  using Cut = typename StructDependencies::Cut;
  using DatabaseEntry = typename Database::database_entry_t;

  struct StructParams : dependency::default_struct_params {
    static constexpr uint32_t num_vars_sign = Params::kMaxNumLeaves;
    static constexpr uint32_t max_cuts_size = Params::kMaxCutsSize;
    static constexpr uint32_t max_reconv = Params::kMaxCutReconv;
  };

  struct ResynCandidate {
    math::static_vector<Signal, Params::kMaxCutsSize> leaves;
    DatabaseEntry entry;
    double reward;
    double arrival_time;
  };

 public:
  RewriteDelayImpl(Ntk& ntk, Database& database, Params const& ps, Stats& st)
      :  // parameters
        ps_(ps),
        // statistics
        st_(st),
        // attributes
        ntk_(ntk),
        database_(database),
        window_manager_(ntk_, ps, st),
        window_simulator_(ntk_),
        dependencies_(ntk_) {}

  void Run();

 private:
  bool SkipNode(Node) const;
  bool TryRewrite(Node);

  double Area(Cut&);
  double RecursiveRef(Node);
  double RecursiveDeref(Node);

 private:
  // Parameters
  Params const& ps_;
  // Statistics
  Stats& st_;
  // Attributes
  Ntk& ntk_;
  WindowManagerT window_manager_;
  WindowSimulatorT window_simulator_;
  StructDependencies dependencies_;
  Database& database_;
};

// Implementation of Run split into steps
template <typename Ntk, typename Database, typename Params, typename Stats>
void RewriteDelayImpl<Ntk, Database, Params, Stats>::Run() {
  auto& times = window_manager_.GetTimes();
  times.Run();
  int32_t const last_tier = times.GetWorstTier();

  for (int32_t tier = last_tier; tier > 0; --tier) {
    bool update = false;
    times.ForeachNode(tier, [&](auto const& n) {
      if (SkipNode(n)) return true;

      window_manager_.Run(n);
      window_simulator_.Run(window_manager_);

      if (!window_manager_.IsValid()) return true;

      if (TryRewrite(n)) {
        st_.RecordRewrite();
        update = true;
      }
      return true;
    });

    // update the levels if there was an optimization
    if (update) times.Run();
  }
}

template <typename Ntk, typename Database, typename Params, typename Stats>
bool RewriteDelayImpl<Ntk, Database, Params, Stats>::SkipNode(Node n) const {
  if (ntk_.fanout_size(n) <= 0) return true;
  if (ntk_.is_pi(n) || ntk_.is_constant(n)) return true;
  if (ntk_.is_dead(n)) return true;
  return false;
}

template <typename Ntk, typename Database, typename Params, typename Stats>
bool RewriteDelayImpl<Ntk, Database, Params, Stats>::TryRewrite(Node n) {
  if (ntk_.num_outputs(n) > 1) return false;
  std::optional<ResynCandidate> best_candidate;
  auto& times = window_manager_.GetTimes();
  times.Run();
  math::static_vector<double, Params::kMaxCutsSize> input_times;
  math::static_vector<Signal, Params::kMaxCutsSize> input_signals;
  math::static_vector<bool, Params::kMaxCutsSize> valid;
  dependencies_.run(window_manager_, window_simulator_);
  dependencies_.foreach_cut([&](auto& cut, auto) {
    double const area_curr = Area(cut);
    auto tt = cut.func[0]._bits;
    input_times.clear();
    input_signals.clear();
    valid.clear();

    for (auto const& f : cut.leaves) {
      input_times.push_back(times.GetArrival(f));
      input_signals.push_back(f);
      valid.push_back(true);
    }
    while (input_times.size() < Params::kMaxCutsSize) {
      input_times.push_back(std::numeric_limits<double>::max());
      input_signals.push_back(Signal{std::numeric_limits<uint32_t>::max()});
      valid.push_back(false);
    }
    auto match =
        database_.boolean_matching(tt, input_times, input_signals, valid);
    if (match) {
      double reward = 0;
      database_.foreach_entry(*match, [&](auto& entry) {
        double arrival_time = 0.0;
        int j = 0;
        for (int i = 0; i < Params::kMaxCutsSize; ++i) {
          if (!valid[i]) continue;
          double const input_time = input_times[i];
          double const delay = entry.delays[j++];
          arrival_time = std::max(arrival_time, input_time + delay);
        }

        reward = times.GetArrival(cut.root) - arrival_time;

        if (ps_.area_constrained && (entry.area > area_curr)) return;

        if (reward > 0 &&
            (!best_candidate || (reward > best_candidate->reward) ||
             ((reward == best_candidate->reward) && (entry.area < area_curr))))
          best_candidate = ResynCandidate{input_signals, entry, reward};
      });
    }
  });

  if (best_candidate) {
    auto const nnew =
        database_.write(best_candidate->entry, ntk_, best_candidate->leaves);
    ntk_.substitute_node(n, ntk_.make_signal(nnew));
    return true;
  }
  return false;
}

template <typename Ntk, typename Database, typename Params, typename Stats>
double RewriteDelayImpl<Ntk, Database, Params, Stats>::Area(Cut& cut) {
  for (auto const& f : cut.leaves) {
    ntk_.incr_fanout_size(ntk_.get_node(f));
  }
  double const cost_curr = RecursiveDeref(cut.root);
  double const cost_test = RecursiveRef(cut.root);
  for (auto const& f : cut.leaves) {
    ntk_.decr_fanout_size(ntk_.get_node(f));
  }

  return cost_curr;
}

template <typename Ntk, typename Database, typename Params, typename Stats>
double RewriteDelayImpl<Ntk, Database, Params, Stats>::RecursiveDeref(Node n) {
  /* terminate? */
  if (ntk_.is_constant(n) || ntk_.is_pi(n)) return 0.0;

  /* recursively collect nodes */
  double area = ntk_.get_area(n);
  ntk_.foreach_fanin(n, [&](auto const& fi) {
    Node const ni = ntk_.get_node(fi);
    if (ntk_.decr_fanout_size(ni) == 0) {
      area += RecursiveDeref(ni);
    }
  });
  return area;
}

template <typename Ntk, typename Database, typename Params, typename Stats>
double RewriteDelayImpl<Ntk, Database, Params, Stats>::RecursiveRef(Node n) {
  /* terminate? */
  if (ntk_.is_constant(n) || ntk_.is_pi(n)) return 0.0;

  /* recursively collect nodes */
  double area = ntk_.get_area(n);

  ntk_.foreach_fanin(n, [&](auto const& fi) {
    Node const ni = ntk_.get_node(fi);
    if (ntk_.incr_fanout_size(ni) == 0) {
      area += RecursiveRef(ni);
    }
  });
  return area;
}

}  // namespace detail

template <class Ntk, typename Db, typename Params,
          typename Stats = RewriteDelayStatsDefault>
void RewriteDelay(Ntk& ntk, Db& db, Params const& ps, Stats* pst = nullptr) {
  Stats st;
  detail::RewriteDelayImpl<Ntk, Db, Params, Stats> p(ntk, db, ps, st);
  p.Run();
  if (pst) *pst = st;
}

}  // namespace cell
}  // namespace opto
}  // namespace rinox
#endif