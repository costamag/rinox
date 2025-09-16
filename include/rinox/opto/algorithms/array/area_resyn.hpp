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

#ifndef RINOX_ALGORITHMS_ARRAY_AREA_RESYN_HPP_
#define RINOX_ALGORITHMS_ARRAY_AREA_RESYN_HPP_

#include "../../../network/visualize.hpp"
#include "../../../windowing/window_manager.hpp"
#include "../../../windowing/window_simulator.hpp"
#include "../../../boolean/spfd.hpp"
#include "../../../boolean/simd.hpp"
#include "../../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../../analyzers/trackers/required_times_tracker.hpp"
#include "../../profilers/area_profiler.hpp"
#include "../../profilers/delay_profiler.hpp"
#include "../../profilers/power_profiler.hpp"
#include "../../profilers/profilers_utils.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <chrono>
#include <optional>

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
namespace algorithms {
namespace array {

struct AreaResynStats {
  mockturtle::stopwatch<>::duration time_window_constr{0};
  mockturtle::stopwatch<>::duration time_window_simula{0};
  mockturtle::stopwatch<>::duration time_0_resub{0};

  double GetTimeWindowConstr() const {  return mockturtle::to_seconds( time_window_constr ); }
  double GetTimeWindowSimula() const {  return mockturtle::to_seconds( time_window_simula ); }
  double GetTime0Resub() const {  return mockturtle::to_seconds( time_0_resub ); }
  double GetTimeTotal() const { return  GetTimeWindowConstr() +
                                        GetTimeWindowSimula() +
                                        GetTime0Resub(); }

  void Report() const {
    std::cout << fmt::format( "    time window construction    = {:5.2f}\n", GetTimeWindowConstr() );
    std::cout << fmt::format( "    time window simulation      = {:5.2f}\n", GetTimeWindowSimula() );
    std::cout << fmt::format( "    time 0 resubstitution       = {:5.2f}\n", GetTime0Resub() );
    std::cout << fmt::format( "    time total                  = {:5.2f}\n", GetTimeTotal() );
  }
};

struct DefaultAreaResynParams {
  // static parameters
  static constexpr uint32_t kMaxNumLeaves = RINOX_MAX_NUM_LEAVES;
  static constexpr uint32_t kNumVarsSign = RINOX_NUM_VARS_SIGN;
  static constexpr uint32_t kMaxCutsSize = RINOX_MAX_CUTS_SIZE;
  static constexpr uint32_t kMaxCubeSpfd = RINOX_MAX_CUBE_SPFD;
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr uint32_t kMaxInpsSolo = 6u;

  // dynamic parameters (mutable, no k-prefix)
  uint32_t stop_after_nth = std::numeric_limits<uint32_t>::max();
  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 10;
  int32_t odc_levels = 0;
  bool preserve_times = true;
  bool use_required = true;
};

namespace detail {

template <class Ntk, typename AreaResynParams = DefaultAreaResynParams>
class AreaResynImpl {
  struct WindowManagerParams {
    static constexpr uint32_t max_num_leaves =
        AreaResynParams::kMaxNumLeaves;
    uint32_t max_num_divisors;
    uint32_t max_reconvergence;
    int32_t odc_levels;
  };

  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using Clock = std::chrono::high_resolution_clock;
  using ArrivalTracker = analyzers::trackers::arrival_times_tracker<Ntk>;
  using RequiredTracker = analyzers::trackers::required_times_tracker<Ntk>;
  using WindowManager =
      windowing::window_manager<Ntk, WindowManagerParams>;
  using WindowManagerStats = windowing::window_manager_stats;
  using WindowSimulator =
      windowing::window_simulator<Ntk, AreaResynParams::kMaxNumLeaves>;

  struct SharedOptoInfo {
    std::unique_ptr<RequiredTracker> required;
    std::unique_ptr<ArrivalTracker> arrivals;
    std::unique_ptr<WindowManager> manager;
    std::unique_ptr<WindowSimulator> simulator;
    Node root;
  };
  struct ReturnOptoData {
    Node old_node;
    Signal new_signal;
  };
 public:
  AreaResynImpl(Ntk& ntk, AreaResynParams ps, AreaResynStats& st)
      : ps_(ps), st_(st), ntk_(ntk) {}

void Run() {
  // ======================================================= //
  // Phase 0: Setup optimization helpers
  // ======================================================= //
  // Initialize arrival time tracker for ntk_. Required for
  // evaluating timing constraints during window construction.
  if ( ps_.preserve_times )
  {
    info_.arrivals = std::make_unique<ArrivalTracker>(ntk_);
    double const worst_delay = info_.arrivals->worst_delay();
    if ( ps_.use_required )
      info_.required = std::make_unique<RequiredTracker>(ntk_, worst_delay);
  }


  // Window manager configuration:
  // - max_num_divisors: upper bound on divisor candidates considered
  // - max_reconvergence: cap on reconvergent paths allowed
  // - odc_levels: depth for observing don't-cares
  WindowManagerParams window_manager_ps{
      ps_.max_num_divisors,
      ps_.max_reconvergence,
      ps_.odc_levels};

  WindowManagerStats window_manager_st;
  info_.manager = std::make_unique<WindowManager>( ntk_, window_manager_ps, window_manager_st );
  // Simulator bound to windows generated by WindowManager.
  // Evaluates candidate replacements for functional equivalence.
  info_.simulator = std::make_unique<WindowSimulator>( ntk_ );

  // ======================================================= //
  // Phase 1: Iterate network gates
  // ======================================================= //
  // For each logic gate (excluding constants and PIs),
  // attempt to resynthesize it by analyzing its fanin window.
  auto const size_init = ntk_.num_gates();
  ntk_.foreach_gate([&](auto const& n) {

    info_.root = n;

    // Skip trivial/unresynthesizable nodes.
    if ( n > size_init || ntk_.is_dead(n) || ntk_.is_pi(n)) return true;
    if ( ntk_.fanin_size( n ) == 1 ) return true;
    // ====================================================== //
    // Phase 2: Window analysis
    // ====================================================== //
    // Construct a localized subgraph (window) around node n.
    // Stopwatch captures cumulative time across all nodes.
    mockturtle::call_with_stopwatch(st_.time_window_constr,
                                    [&]() { info_.manager->run(n); });

    // ====================================================== //
    // Phase 3: Simulation of candidate functions
    // ====================================================== //
    // Evaluate functionality of the window using fast bit-parallel
    // simulation. This detects equivalences and filters candidates
    // before invoking cost-driven resubstitution.
    mockturtle::call_with_stopwatch(st_.time_window_simula,
                                    [&]() { info_.simulator->run(*info_.manager); });

    // ====================================================== //
    // Phase 4: Resubstitution (not yet implemented)
    // ====================================================== //
    // TODO(andrea): integrate divisor enumeration and cost-driven
    // replacement strategy. Initial phase may attempt 0-resub
    // (reuse of existing functionally equivalent nodes).
    auto resub0 = mockturtle::call_with_stopwatch(st_.time_0_resub,
                              [&]() { return Try0Resub(); });
    if ( resub0 ) { ApplyOpto( *resub0 ); return true; }
    return true;  // continue traversal
  });
}

private:
  void ApplyOpto( ReturnOptoData const& opto )
  {
    ntk_.substitute_node( opto.old_node, opto.new_signal );
  }

  std::optional<ReturnOptoData> Try0Resub()
  {
    if ( ntk_.num_outputs( info_.root ) != 1 )
      return std::nullopt;

    Signal froot = ntk_.make_signal( info_.root );
    auto const& signature_root = info_.simulator->Get( froot );

    std::optional<ReturnOptoData> res;

    info_.manager->foreach_divisor( [&]( auto const& f, auto i ){
      if ( res ) return;
      if ( ps_.preserve_times )
      {
        double const time_max = ps_.use_required ? info_.required->get_time( froot ) : info_.arrivals->get_time( froot );
        if ( info_.arrivals->get_time( froot ) > time_max ) return;
      }
      auto const& signature_div = info_.simulator->Get( f );
      auto const diff = boolean::binary_xor( signature_div, signature_root );

      if ( ps_.odc_levels <= 0 )
      {
        auto const num_ones = kitty::count_ones( diff );
        if ( num_ones == 0 )
          res = { info_.root, f };
        else if ( num_ones == diff.num_bits() )
          res = { info_.root, ntk_.create_not(f) };
      }
      else
      {
        auto const& care = info_.simulator->get_careset();
        auto const num_ones = kitty::count_ones( boolean::binary_and( diff, care ) );
        if ( num_ones == 0 )
          res = { info_.root, f };
        else if ( num_ones == kitty::count_ones( care ) )
          res = { info_.root, ntk_.create_not(f) };
      }
    } );
    return res;
  }

 private:
  Ntk& ntk_;
  AreaResynParams ps_;
  AreaResynStats& st_;
  SharedOptoInfo info_;
};

}  // namespace detail

template <class Ntk, typename AreaResynParams = DefaultAreaResynParams>
void AreaResyn(Ntk& ntk, AreaResynParams ps = {},
               AreaResynStats* pst = nullptr) {
  AreaResynStats st;
  detail::AreaResynImpl<Ntk, AreaResynParams> p(ntk, ps, st);

  p.Run();
  if (pst != nullptr) *pst = st;
}

}  // namespace array
}  // namespace algorithms
}  // namespace opto
}  // namespace rinox

#endif  // RINOX_ALGORITHMS_ARRAY_AREA_RESYN_HPP_
