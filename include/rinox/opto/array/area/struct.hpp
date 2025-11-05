/* rinox: C++ logic network library
 * Copyright (C) 2025 EPFL
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
  \file struct.hpp
  \brief Structural resinthesis of LUT networks

  \author Andrea Costamagna
*/

#ifndef RINOX_OPTO_ARRAY_AREA_STRUCT_HPP_
#define RINOX_OPTO_ARRAY_AREA_STRUCT_HPP_

#include <fmt/format.h>

#include <chrono>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <optional>

#include "../../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../../boolean/spfd.hpp"
#include "../../../windowing/window_manager.hpp"
#include "../../../windowing/window_simulator.hpp"
#include "../../../synthesis/lut/synth.hpp"

namespace rinox {
namespace opto {
namespace array {

struct AreaStructStaticParams {
  static constexpr uint32_t kMaxNumLeaves = 15u;
  static constexpr uint32_t kMaxInpsSolo = 6u;
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr synthesis::lut::LutSynthMethod kLutSynth =
    synthesis::lut::LutSynthMethod::kSpinDecSolo;
};

struct AreaStructDynamicParams {
  uint32_t stop_after_nth = std::numeric_limits<uint32_t>::max();
  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 10;
  int32_t odc_levels = 2;
};

struct AreaStructStats {
  mockturtle::stopwatch<>::duration time_window_constr{0};
  mockturtle::stopwatch<>::duration time_window_simula{0};
  mockturtle::stopwatch<>::duration time_resyn_solo{0};
  mockturtle::stopwatch<>::duration time_resyn_dual{0};
  uint32_t num_resyn_solo{0};
  uint32_t num_resyn_dual{0};

  double GetTimeWindowConstr() const { return mockturtle::to_seconds(time_window_constr); }
  double GetTimeWindowSimula() const { return mockturtle::to_seconds(time_window_simula); }
  double GetTimeTotal() const {
      return  GetTimeWindowConstr() +
              GetTimeWindowSimula();
  }

  void RecordResynSolo() { num_resyn_solo++; }
  void RecordResynDual() { num_resyn_dual++; }

  uint32_t GetNumResynSolo() const { return num_resyn_solo; }
  uint32_t GetNumResynDual() const { return num_resyn_dual; }
  uint32_t GetNumResynTotal() const { return  GetNumResynSolo() +
                                            GetNumResynDual(); }

  void Report() const {
    std::cout << fmt::format("================ TIME =====================\n");
    std::cout << fmt::format("    time window construction      = {:5.2f}\n",
                             GetTimeWindowConstr());
    std::cout << fmt::format("    time window simulation        = {:5.2f}\n",
                             GetTimeWindowSimula());

    std::cout << fmt::format("    time total                    = {:5.2f}\n",
                             GetTimeTotal());
    std::cout << fmt::format("================ COUNT ====================\n");
    std::cout << fmt::format("    num one-output resynthesis   = {:5d}\n",
                             GetNumResynSolo());
    std::cout << fmt::format("    num two-output resynthesis   = {:5d}\n",
                             GetNumResynDual());

    std::cout << fmt::format("    num total                    = {:5d}\n",
                             GetNumResynTotal());
  }
};

namespace detail {

template <class Ntk, typename StaticParams = AreaStructStaticParams>
class AreaStructImpl {
  struct WindowManagerParams {
    static constexpr uint32_t max_num_leaves = StaticParams::kMaxNumLeaves;
    uint32_t max_num_divisors;
    uint32_t max_reconvergence;
    int32_t odc_levels;
  };

  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using ArrivalsTracker = analyzers::trackers::arrival_times_tracker<Ntk>;
  using Signature = kitty::static_truth_table<StaticParams::kMaxNumLeaves>;
  using WindowManager = windowing::window_manager<Ntk, WindowManagerParams>;
  using WindowSimulator = windowing::window_simulator<Ntk, StaticParams::kMaxNumLeaves>;
  static constexpr uint32_t kMaxNumOutputs = Ntk::kMaxNumOutputs;
  using Spfds = boolean::spfds_manager<Signature, 1u << StaticParams::kMaxInpsSolo, kMaxNumOutputs>;

  struct LutSynthParams : synthesis::lut::DefaultLutSynthParams {
    static constexpr uint32_t kNumVarsCube = StaticParams::kMaxNumLeaves;
    static constexpr uint32_t kMaxInpsSolo = StaticParams::kMaxInpsSolo;
    static constexpr uint32_t kMaxInpsDual = StaticParams::kMaxInpsDual;
    static constexpr uint32_t kMaxNumOutputs = Ntk::kMaxNumOutputs;
    static constexpr uint32_t kMaxNumGates = 32u;
  };

  using Synth = synthesis::lut::LutSynthDispatch<synthesis::lut::LutSynthMethod::kSpinDecSolo, Ntk, WindowManager, WindowSimulator, LutSynthParams>;

 public:
  AreaStructImpl(Ntk& ntk, AreaStructDynamicParams ps, AreaStructStats& st)
      : ntk_(ntk),
        ps_(ps),
        st_(st),
        arrivals_(ntk_) 
  {
    for (int n = 0; n <= StaticParams::kMaxInpsSolo; ++n) {
      proj_fns_.emplace_back();
      for (int i = 0; i < n; ++i) {
        proj_fns_[n].emplace_back(n);
        kitty::create_nth_var(proj_fns_[n][i], i);
      }
    }
  }

  void Run();
  void AnalyzeStructure();

 private:
  bool SkipNode(Node) const;
  bool ResynSolo(Node, WindowManager& , WindowSimulator&);

 private:
  Ntk& ntk_;
  AreaStructDynamicParams ps_;
  AreaStructStats& st_;

  Spfds spfd_;
  Synth synth_;
  ArrivalsTracker arrivals_;
  std::vector<std::vector<kitty::dynamic_truth_table>> proj_fns_;

};

// Implementation of Run split into steps
template <class Ntk, typename StaticParams>
void AreaStructImpl<Ntk, StaticParams>::Run() {
  const auto nmax = ntk_.size();
  int cnt = 0;
  AnalyzeStructure();
  // windwong infomation
  WindowManagerParams ps_manager;
  windowing::window_manager_stats st_manager;
  ps_manager.max_num_divisors = ps_.max_num_divisors;
  ps_manager.max_reconvergence = ps_.max_reconvergence;
  ps_manager.odc_levels = ps_.odc_levels;

  WindowManager manager(ntk_, ps_manager, st_manager);
  WindowSimulator simulator(ntk_);

  ntk_.foreach_gate([&](auto const& n) {
    if (ntk_.is_pi(n)) return true;
    if (st_.GetNumResynTotal() >= ps_.stop_after_nth) return false;
    if (SkipNode(n)) return true;
    if (ntk_.num_outputs(n) != 1) return true;

    mockturtle::call_with_stopwatch(st_.time_window_constr, [&]() { manager.run(n); });
    mockturtle::call_with_stopwatch(st_.time_window_simula, [&]() { simulator.run(manager);});

    if ( mockturtle::call_with_stopwatch(st_.time_resyn_solo, [&]() { return ResynSolo(n, manager, simulator); }))
    {
      st_.RecordResynSolo();
    }
    return true;
  });
}

template <class Ntk, typename StaticParams>
void AreaStructImpl<Ntk, StaticParams>::AnalyzeStructure() {

  // row leaves column mffc size
  std::vector<std::vector<uint32_t>> distr;

  const auto nmax = ntk_.size();
  int cnt = 0;

  // windwong infomation
  WindowManagerParams ps_manager;
  windowing::window_manager_stats st_manager;
  ps_manager.max_num_divisors = ps_.max_num_divisors;
  ps_manager.max_reconvergence = ps_.max_reconvergence;
  ps_manager.odc_levels = ps_.odc_levels;

  WindowManager manager(ntk_, ps_manager, st_manager);
  WindowSimulator simulator(ntk_);

  ntk_.foreach_gate([&](auto const& n) {
    if (ntk_.is_pi(n)) return true;
    if (SkipNode(n)) return true;
    if (ntk_.num_outputs(n) != 1) return true;

    manager.run(n);
    simulator.run(manager);

    while( distr.size() <= manager.num_inputs() )
      distr.emplace_back();
    while( distr[ manager.num_inputs() ].size() <= manager.num_mffc() )
      distr[ manager.num_inputs() ].push_back(0);
    distr[ manager.num_inputs() ][manager.num_mffc()]++;

    return true;
  });

  for ( uint32_t i = 1; i < distr.size(); ++i )
  {
    for ( uint32_t j = 2; j < distr[i].size(); ++j )
    {
      std::cout << "(" << i << ", " << j << "):" << distr[i][j] << " ";
    }
    std::cout << std::endl;
  }

}

template <class Ntk, typename StaticParams>
bool AreaStructImpl<Ntk, StaticParams>::SkipNode(Node n) const {
  if (ntk_.fanout_size(n) <= 0) return true;
  if (ntk_.is_pi(n) || ntk_.is_constant(n)) return true;
  if (ntk_.is_dead(n)) return true;
  return false;
}

template <class Ntk, typename StaticParams>
bool AreaStructImpl<Ntk, StaticParams>::ResynSolo(Node n, WindowManager& manager, WindowSimulator& simulator ) {
  uint32_t const K = manager.num_inputs();
  uint32_t const M = manager.num_mffc();
  uint32_t const T = manager.num_tfos();
  uint32_t const O = manager.num_outputs();
  static constexpr uint32_t k = StaticParams::kMaxInpsSolo;
  Signature const& care = simulator.get_careset();
  
  uint32_t num_dcs = (1u << StaticParams::kMaxNumLeaves)-kitty::count_ones(care);
  std::cout << fmt::format("Processing node {} with {} inputs {} MFFC {} TFOs {} Outouts and {} DCs\n", n, K, M, T, O, num_dcs);
  if ( (num_dcs == 0 ) && ( (K==0) || M <= 1 || (K > ( (M-1)*k - (M-2) )) )) return false;
  
  //bool const success = synth_.Run(manager, simulator);
  //if(success)
  //{
  //  auto const& sub_problems = synth_.GetSubproblems();
  //  auto const fnew = Insert(sub_problems);
  //  ntk_.substitute_node(n, *fnew);
  //  return true;
  //}
  return false;
}

}  // namespace detail

// Public API
template <class Ntk, typename StaticParams = AreaStructStaticParams>
void AreaStruct(Ntk& ntk, AreaStructDynamicParams ps = {},
               AreaStructStats* pst = nullptr) {
  AreaStructStats st;
  detail::AreaStructImpl<Ntk, StaticParams> p(ntk, ps, st);
  p.Run();
  if (pst) *pst = st;
}

}  // namespace array
}  // namespace opto
}  // namespace rinox

#endif  // RINOX_OPTO_ARRAY_AREA_STRUCT_HPP_
