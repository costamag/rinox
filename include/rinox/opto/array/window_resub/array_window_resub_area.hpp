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

#ifndef RINOX_ALGORITHMS_ARRAY_WINDOW_RESUB_AREA_HPP_
#define RINOX_ALGORITHMS_ARRAY_WINDOW_RESUB_AREA_HPP_

#include <fmt/format.h>

#include <chrono>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <optional>

#include "../../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../../analyzers/evaluators/levels_evaluator.hpp"
#include "../../../boolean/simd.hpp"
#include "../../../boolean/spfd.hpp"
#include "../../../dependency/window/oracles.hpp"
#include "../../../network/visualize.hpp"
#include "../../../synthesis/lut/synth.hpp"
#include "../../../windowing/simula_manager.hpp"
#include "../../../windowing/window_manager.hpp"
#include "../../../windowing/window_simulator.hpp"
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
namespace array {

struct WindowResubAreaParamsDefault {
  static constexpr uint32_t kMaxNumLeaves = RINOX_MAX_NUM_LEAVES;
  static constexpr uint32_t kMaxCutsSize = RINOX_MAX_CUTS_SIZE;
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr uint32_t kMaxInpsSolo = 6u;
  uint32_t stop_after_nth = std::numeric_limits<uint32_t>::max();
  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 10;
  int32_t odc_levels = 0;
  bool try_0resub = false;
  bool try_1resub = false;
  bool try_2resub = false;
};

struct WindowResubAreaStatsDefault {
  mockturtle::stopwatch<>::duration time_window_construction{0};
  mockturtle::stopwatch<>::duration time_window_simulation{0};
  mockturtle::stopwatch<>::duration time_0resub{0};
  mockturtle::stopwatch<>::duration time_1resub{0};
  mockturtle::stopwatch<>::duration time_2resub{0};
  uint32_t num_changes{0};
  uint32_t num_0resub{0};
  uint32_t num_1resub{0};
  uint32_t num_2resub{0};

  double GetTimeWindowConstr() const { return mockturtle::to_seconds(time_window_construction); }
  double GetTimeWindowSimula() const { return mockturtle::to_seconds(time_window_simulation); }
  double GetTime0Resub() const {return mockturtle::to_seconds(time_0resub);}
  double GetTime1Resub() const {return mockturtle::to_seconds(time_1resub);}
  double GetTime2Resub() const {return mockturtle::to_seconds(time_2resub);}

  double GetTimeTotal() const { return  GetTimeWindowConstr() +
                                        GetTimeWindowSimula() +
                                        GetTimeOpto(); }

  double GetTimeOpto() const { return GetTime0Resub() +
                                      GetTime1Resub() +
                                      GetTime2Resub(); }

  void Record0Resub() { num_0resub++; }
  void Record1Resub() { num_1resub++; }
  void Record2Resub() { num_2resub++; }

  uint32_t GetNum0Resub() const { return num_0resub; }
  uint32_t GetNum1Resub() const { return num_1resub; }
  uint32_t GetNum2Resub() const { return num_2resub; }

  uint32_t GetNumOptos() const { return GetNum0Resub() +
                                        GetNum1Resub() +
                                        GetNum2Resub(); }

  void Report() const {
    std::cout << fmt::format("================ TIME  ====================\n");
    std::cout << fmt::format("    time window construction      = {:5.2f}\n",
                             GetTimeWindowConstr());
    std::cout << fmt::format("    time window simulation        = {:5.2f}\n",
                             GetTimeWindowSimula());

    std::cout << fmt::format("    time 0-resub                  = {:5.2f}\n",
                             GetTime0Resub());
    std::cout << fmt::format("    time 1-resub                  = {:5.2f}\n",
                             GetTime1Resub());
    std::cout << fmt::format("    time 2-resub                  = {:5.2f}\n",
                             GetTime2Resub());

    std::cout << fmt::format("    time total                    = {:5.2f}\n",
                             GetTimeTotal());
    std::cout << fmt::format("================ QUALITY ==================\n");
    std::cout << fmt::format("    count 0-resub                 = {:5d}\n",
                             GetNum0Resub());
    std::cout << fmt::format("    count 1-resub                 = {:5d}\n",
                             GetNum1Resub());
    std::cout << fmt::format("    count 2-resub                 = {:5d}\n",
                             GetNum2Resub());
    std::cout << fmt::format("    count total                   = {:5d}\n",
                             GetNumOptos());
  }
};

template <typename Ntk>
struct ReturnOptoData {
  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;

  Node old_node;
  Signal new_signal;
};

namespace detail {

template <class Ntk, typename Params, typename Stats>
class WindowResubAreaImpl {

  struct WindowManagerParams {
    static constexpr uint32_t max_num_leaves = Params::kMaxNumLeaves;
    uint32_t max_num_divisors;
    uint32_t max_reconvergence;
    int32_t odc_levels;
  };

  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using Signature = kitty::static_truth_table<Params::kMaxNumLeaves>;
  using WindowManager = windowing::window_manager<Ntk, WindowManagerParams>;
  using WindowSimulator = windowing::window_simulator<Ntk, Params::kMaxNumLeaves>;
  using Spfds = boolean::spfd<Signature, 1u << Params::kMaxCutsSize>;
  static constexpr uint32_t kMaxNumOutputs = Ntk::kMaxNumOutputs;

 public:
  WindowResubAreaImpl(Ntk& ntk, Params ps, Stats& st)
      : // parameters
        ps_(ps),
        window_manager_ps_{ps_.max_num_divisors, ps_.max_reconvergence, ps_.odc_levels},
        // statistics
        st_(st),
        // attributes
        ntk_(ntk),
        window_manager_(ntk_, window_manager_ps_, window_manager_st_),
        window_simulator_(ntk_) {
  }

  void Report();
  void Run();

 private:
  bool SkipNode(Node) const;
  bool Try0Resub(Node);
  bool Try1Resub(Node);
  bool Try2Resub(Node);

 private:
  // Parameters
  Params ps_;
  WindowManagerParams window_manager_ps_;
  // Statistics
  Stats& st_;
  windowing::window_manager_stats window_manager_st_;
  windowing::SimulaManagerStats window_simulator_st_;
 // Attributes
  Ntk& ntk_;
  Spfds spfd_;
  WindowManager window_manager_;
  WindowSimulator window_simulator_;
};

template <typename Ntk, typename Params, typename Stats>
void WindowResubAreaImpl<Ntk, Params, Stats>::Report() {
  std::vector<uint32_t> mffc_sizes_solo;
  std::vector<uint32_t> mffc_sizes_dual;
  ntk_.foreach_gate([&](auto const& n) {
    window_manager_.run(n);
    uint32_t mffc_size = window_manager_.num_mffc();
    if (ntk_.num_outputs(n)==1)
    {
      while (mffc_sizes_solo.size() <= mffc_size) mffc_sizes_solo.push_back(0);
      mffc_sizes_solo[mffc_size] += 1;
    }
    else
    {
      while (mffc_sizes_dual.size() <= mffc_size) mffc_sizes_dual.push_back(0);
      mffc_sizes_dual[mffc_size] += 1;
    }
    return true;
  });
  for (uint32_t s = 0; s < mffc_sizes_solo.size(); ++s) {
    std::cout << s << " 1 : " << mffc_sizes_solo[s] << std::endl;
  }
  for (uint32_t s = 0; s < mffc_sizes_dual.size(); ++s) {
    std::cout << s << " 2 : " << mffc_sizes_dual[s] << std::endl;
  }
}

// Implementation of Run split into steps
template <typename Ntk, typename Params, typename Stats>
void WindowResubAreaImpl<Ntk, Params, Stats>::Run() {

  rinox::analyzers::evaluators::LevelsEvaluator<Ntk> levels( ntk_ );

  int const last_level = levels.GetLast();

  for(int level=last_level; level>0; --level)
  {
    bool update=false;
    levels.ForeachNode(level, [&](auto const& n){
      if (SkipNode(n)) return true;
      mockturtle::call_with_stopwatch(st_.time_window_construction, [&]() { window_manager_.run(n); });
      mockturtle::call_with_stopwatch(st_.time_window_simulation, [&]() { window_simulator_.run(window_manager_);});
      
      if (ps_.try_0resub && mockturtle::call_with_stopwatch(st_.time_0resub,
                                          [&]() { return Try0Resub(n); })) {
        st_.Record0Resub();
        update=true;
        return true;
      }
      if (ps_.try_1resub && mockturtle::call_with_stopwatch(st_.time_1resub,
                                          [&]() { return Try1Resub(n); })) {
        st_.Record1Resub();
        update=true;
        return true;
      }
      if (ps_.try_2resub && mockturtle::call_with_stopwatch(st_.time_2resub,
                                          [&]() { return Try2Resub(n); })) {
        st_.Record2Resub();
        update=true;
        return true;
      }
      return true;
    });
    // update the levels if there was an optimization
    if(update)
      levels.ComputeLevels();
  }
}

template <typename Ntk, typename Params, typename Stats>
bool WindowResubAreaImpl<Ntk, Params, Stats>::SkipNode(Node n) const {
  if (ntk_.fanout_size(n) <= 0) return true;
  if (ntk_.is_pi(n) || ntk_.is_constant(n)) return true;
  if (ntk_.is_dead(n)) return true;
  return false;
}

template <typename Ntk, typename Params, typename Stats>
bool WindowResubAreaImpl<Ntk, Params, Stats>::Try0Resub(Node n) {
  std::vector<Signal> new_signals;

  ntk_.foreach_output(n, [&](auto const& f){
    auto const& onset = window_simulator_.Get(f);
    std::optional<Signal> resub;
    window_manager_.foreach_divisor([&](auto const& d, auto) {
      if(!resub)
      {
        auto const& divisor = window_simulator_.Get(d);
        if (kitty::equal(onset, divisor))
        {
          kitty::print_binary(divisor);
          std::cout << " " << d << std::endl;
          resub = d;
        }
      }
    });
    if (resub) new_signals.push_back(*resub);
  });
  if(new_signals.size()==ntk_.num_outputs(n))
  {
    ntk_.substitute_node(n, new_signals);
    return true;
  }
  return false;
}

template <typename Ntk, typename Params, typename Stats>
bool WindowResubAreaImpl<Ntk, Params, Stats>::Try1Resub(Node n) {
  return false;
}

template <typename Ntk, typename Params, typename Stats>
bool WindowResubAreaImpl<Ntk, Params, Stats>::Try2Resub(Node n) {
  return false;
}

}  // namespace detail

// Public API
template <class Ntk, typename Params = WindowResubAreaParamsDefault, typename Stats = WindowResubAreaStatsDefault>
void WindowResubArea( Ntk& ntk,
                      Params ps = {},
                      Stats* pst = nullptr) {
  Stats st;
  detail::WindowResubAreaImpl<Ntk, Params, Stats> p(ntk, ps, st);
  p.Run();
  st.Report();
  if (pst) *pst = st;
}

}  // namespace array
}  // namespace opto
}  // namespace rinox

#endif  // RINOX_ALGORITHMS_ARRAY_AREA_ORCHESTRA_HPP_
