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

#ifndef RINOX_ALGORITHMS_CELL_REWRITE_AREA_HPP_
#define RINOX_ALGORITHMS_CELL_REWRITE_AREA_HPP_

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
#include "../../../dependency/struct_dependencies.hpp"
#include "../../../boolean/simd.hpp"
#include "../../../boolean/spfd.hpp"
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

struct RewriteAreaParamsDefault : windowing::WindowManagerParamsDefault {
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
  int32_t max_num_tiers_tfi = 10;
  int32_t max_num_tiers_tfo = 0;
};

struct RewriteAreaStatsDefault : windowing::WindowManagerStatsDefault {
  mockturtle::stopwatch<>::duration time_rewrite{0};
  uint32_t num_rewrite{0};

  double GetTimeRewrite() const { return mockturtle::to_seconds(time_rewrite); }

  double GetTimeTotal() const {
    return GetTimeRewrite();
  }

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
class RewriteAreaImpl {

  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using Signature = kitty::static_truth_table<Params::kMaxNumLeaves>;
  using WindowManagerT = windowing::WindowManager<Ntk, Params>;
  using WindowSimulatorT =
      windowing::WindowSimulator<Ntk, Params::kMaxNumLeaves>;
  static constexpr uint32_t kMaxNumOutputs = Ntk::kMaxNumOutputs;
  using TimesEvaluatorT = typename analyzers::evaluators::TimesEvaluator<Ntk>;

  struct StructParams : dependency::default_struct_params
  {
    static constexpr uint32_t num_vars_sign = Params::kMaxNumLeaves;
    static constexpr uint32_t max_cuts_size = Params::kMaxCutsSize;
    static constexpr uint32_t max_reconv = Params::kMaxCutReconv;
  };

  using StructDependencies = dependency::struct_dependencies<Ntk, StructParams>;

 public:
  RewriteAreaImpl(Ntk& ntk, Database & database, Params ps, Stats& st)
      :  // parameters
        ps_(ps),
        // statistics
        st_(st),
        // attributes
        ntk_(ntk),
        database_(database),
        window_manager_(ntk_, ps, st),
        window_simulator_(ntk_),
        dependencies_(ntk_)
        {}

  void Run();

 private:
  bool SkipNode(Node) const;
  bool TryRewrite(Node);

 private:
  // Parameters
  Params ps_;
  // Statistics
  Stats& st_;
  // Attributes
  Ntk& ntk_;
  WindowManagerT window_manager_;
  WindowSimulatorT window_simulator_;
  StructDependencies dependencies_;
  Database & database_;
};

// Implementation of Run split into steps
template <typename Ntk, typename Database, typename Params, typename Stats>
void RewriteAreaImpl<Ntk, Database, Params, Stats>::Run() {
  auto & times = window_manager_.GetTimes();
  times.Run();
  int32_t const last_tier = times.GetWorstTier();

  for (int32_t tier = last_tier; tier > 0; --tier) {
    bool update = false;
    times.ForeachNode(tier, [&](auto const& n) {
      if (SkipNode(n)) return true;

      window_manager_.Run(n);
      window_simulator_.Run(window_manager_);

      if (!window_manager_.IsValid()) return true;

      if( TryRewrite(n) ) {
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
bool RewriteAreaImpl<Ntk, Database, Params, Stats>::SkipNode(Node n) const {
  if (ntk_.fanout_size(n) <= 0) return true;
  if (ntk_.is_pi(n) || ntk_.is_constant(n)) return true;
  if (ntk_.is_dead(n)) return true;
  return false;
}

template <typename Ntk, typename Database, typename Params, typename Stats>
bool RewriteAreaImpl<Ntk, Database, Params, Stats>::TryRewrite(Node n) {

  std::cout << "TFI" << std::endl;
  window_manager_.ForeachSignalDivs([&](auto const& d, auto i){
    std::cout << d.index << " ";
  });
  std::cout << std::endl;

  double best_reward = 0;
  dependencies_.run( window_manager_, window_simulator_ );
  dependencies_.foreach_cut( [&]( auto& cut, auto i ) {
    for(auto x : cut) std::cout << x << " ";
    for(auto const& tt : cut.func)
    {
      kitty::print_binary(tt);
      std::cout << " " << i << " \n";

      std::vector<Signal> loc_leaves;
      loc_leaves.resize( spec.inputs.size() );
      std::vector<double> loc_times;
      loc_times.resize( spec.inputs.size() );
      std::vector<typename decomposer_t::cut_func_t const*> loc_sims = sim_ptrs;
      std::transform( spec.inputs.begin(), spec.inputs.end(), loc_leaves.begin(), [&]( auto const& lit ) { return signals[lit]; } );
      std::transform( spec.inputs.begin(), spec.inputs.end(), loc_times.begin(), [&]( auto const& lit ) { return times[lit]; } );
      // Perform boolean matching
      auto row = database_->boolean_matching( ctt, loc_times, loc_leaves, loc_sims );

      if ( row )
      {
        auto [index, cost_cand] = evaluate( *row, loc_leaves );
        if ( cost_cand < best_loc_cost )
        {
          best_loc_cost = cost_cand;
          best_database_node = std::make_optional( index );
          best_loc_leaves = loc_leaves;
          best_loc_sims = loc_sims;
        }
      }


    }
  } );

  if ( best_reward > 0 )
  {
    st_.num_rewrite++;
    return true;
  }
  return false;
}

}

template <class Ntk, typename Db, typename Params = RewriteAreaParamsDefault,
          typename Stats = RewriteAreaStatsDefault>
void RewriteArea(Ntk& ntk, Db & db, Params ps = {}, Stats* pst = nullptr) {

  Stats st;
  detail::RewriteAreaImpl<Ntk, Db, Params, Stats> p(ntk, db, ps, st);
  p.Run();
  st.Report();
  if (pst) *pst = st;
}

}
}
}
#endif