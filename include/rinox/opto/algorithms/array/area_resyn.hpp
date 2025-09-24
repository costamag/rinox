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

#include <fmt/format.h>

#include <chrono>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <optional>

#include "../../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../../analyzers/trackers/required_times_tracker.hpp"
#include "../../../boolean/simd.hpp"
#include "../../../boolean/spfd.hpp"
#include "../../../dependency/window/oracles.hpp"
#include "../../../network/visualize.hpp"
#include "../../../synthesis/lut/synth.hpp"
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
namespace algorithms {
namespace array {

struct DefaultAreaResynStaticParams {
  static constexpr uint32_t kMaxNumLeaves = RINOX_MAX_NUM_LEAVES;
  static constexpr uint32_t kNumVarsSign = RINOX_NUM_VARS_SIGN;
  static constexpr uint32_t kMaxCutsSize = RINOX_MAX_CUTS_SIZE;
  static constexpr uint32_t kMaxCubeSpfd = RINOX_MAX_CUBE_SPFD;
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr uint32_t kMaxInpsSolo = 6u;
  static constexpr synthesis::lut::LutSynthMethod kLutSynth =
      synthesis::lut::LutSynthMethod::kShannonSolo;
  static constexpr dependency::window::Oracles kSpfdMeth =
      dependency::window::Oracles::kPivs;
};

struct DefaultAreaResynDynamicParams {
  uint32_t stop_after_nth = std::numeric_limits<uint32_t>::max();
  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 10;
  int32_t odc_levels = 0;
  bool preserve_times = true;
  bool use_required = true;
  bool try_0resub_solo = true;
  bool try_rewire_solo = true;
  bool try_struct_solo = true;
  bool try_window_solo = true;
};

struct AreaResynStats {
  mockturtle::stopwatch<>::duration time_window_constr{0};
  mockturtle::stopwatch<>::duration time_window_simula{0};
  mockturtle::stopwatch<>::duration time_0resub_solo{0};
  mockturtle::stopwatch<>::duration time_rewire_solo{0};
  mockturtle::stopwatch<>::duration time_struct_solo{0};
  mockturtle::stopwatch<>::duration time_window_solo{0};
  uint32_t num_changes{0};
  uint32_t num_0resub_solo{0};
  uint32_t num_rewire_solo{0};
  uint32_t num_struct_solo{0};
  uint32_t num_window_solo{0};
  synthesis::lut::LutSynthStats st_synth;

  double GetTimeWindowConstr() const {
    return mockturtle::to_seconds(time_window_constr);
  }
  double GetTimeWindowSimula() const {
    return mockturtle::to_seconds(time_window_simula);
  }
  double GetTime0ResubSolo() const {
    return mockturtle::to_seconds(time_0resub_solo);
  }
  double GetTimeRewireSolo() const {
    return mockturtle::to_seconds(time_rewire_solo);
  }
  double GetTimeStructSolo() const {
    return mockturtle::to_seconds(time_struct_solo);
  }
  double GetTimeWindowSolo() const {
    return mockturtle::to_seconds(time_window_solo);
  }
  double GetTimeTotal() const {
    return GetTimeWindowConstr() + GetTimeWindowSimula() + GetTime0ResubSolo() +
           GetTimeRewireSolo() + GetTimeStructSolo() + GetTimeWindowSolo();
  }

  void Record0ResubSolo() { num_0resub_solo++; }
  void RecordRewireSolo() { num_rewire_solo++; }
  void RecordStructSolo() { num_struct_solo++; }
  void RecordWindowSolo() { num_window_solo++; }

  uint32_t GetNum0ResubSolo() const { return num_0resub_solo; }
  uint32_t GetNumRewireSolo() const { return num_rewire_solo; }
  uint32_t GetNumStructSolo() const { return num_struct_solo; }
  uint32_t GetNumWindowSolo() const { return num_window_solo; }
  uint32_t GetNumChanges() const {
    return GetNum0ResubSolo() + GetNumRewireSolo() + GetNumStructSolo() +
           GetNumWindowSolo();
  }

  void Report() const {
    std::cout << fmt::format("================ TIME  ====================\n");
    std::cout << fmt::format("    time window construction      = {:5.2f}\n",
                             GetTimeWindowConstr());
    std::cout << fmt::format("    time window simulation        = {:5.2f}\n",
                             GetTimeWindowSimula());
    std::cout << fmt::format("    time single-output 0-resub       = {:5.2f}\n",
                             GetTime0ResubSolo());
    std::cout << fmt::format("    time single-output rewiring      = {:5.2f}\n",
                             GetTimeRewireSolo());
    std::cout << fmt::format("    time single-output rewriting     = {:5.2f}\n",
                             GetTimeStructSolo());
    std::cout << fmt::format("    time single-output window-resub  = {:5.2f}\n",
                             GetTimeStructSolo());
    std::cout << fmt::format("    time total                       = {:5.2f}\n",
                             GetTimeTotal());
    std::cout << fmt::format("================ QUALITY ==================\n");
    std::cout << fmt::format("    count single-output resub0       = {:5d}\n",
                             GetNum0ResubSolo());
    std::cout << fmt::format("    count single-output rewire       = {:5d}\n",
                             GetNumRewireSolo());
    std::cout << fmt::format("    count single-output struct       = {:5d}\n",
                             GetNumStructSolo());
    std::cout << fmt::format("    count single-output window resub = {:5d}\n",
                             GetNumWindowSolo());
    std::cout << fmt::format("    count total                      = {:5d}\n",
                             GetNumChanges());
    st_synth.Report();
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

template <class Ntk, typename Params = DefaultAreaResynDynamicParams,
          typename StaticParams = DefaultAreaResynStaticParams>
class AreaResynImpl {
  struct WindowManagerParams {
    static constexpr uint32_t max_num_leaves = StaticParams::kMaxNumLeaves;
    uint32_t max_num_divisors;
    uint32_t max_reconvergence;
    int32_t odc_levels;
  };

  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using ArrivalsTracker = analyzers::trackers::arrival_times_tracker<Ntk>;
  using RequiredTracker = analyzers::trackers::required_times_tracker<Ntk>;
  using Signature = kitty::static_truth_table<StaticParams::kMaxNumLeaves>;
  using WindowManager = windowing::window_manager<Ntk, WindowManagerParams>;
  using WindowSimulator =
      windowing::window_simulator<Ntk, StaticParams::kMaxNumLeaves>;

  struct SharedOptoWindow {
    WindowManager manager;
    WindowSimulator simulator;
    Node root{};
  };

  struct SharedOptoTracks {
    RequiredTracker* required = nullptr;
    ArrivalsTracker* arrivals = nullptr;
  };

  struct LutSynthParams : synthesis::lut::DefaultLutSynthParams {
    static constexpr uint32_t kNumVarsCube = StaticParams::kMaxNumLeaves;
    static constexpr uint32_t kMaxInpsSolo = StaticParams::kMaxInpsSolo;
    static constexpr uint32_t kMaxInpsDual = StaticParams::kMaxInpsDual;
    static constexpr uint32_t kMaxNumGates = 32u;
  };

  using LutSynth =
      typename synthesis::lut::LutSynthDispatch<StaticParams::kLutSynth, Ntk,
                                                WindowManager, WindowSimulator,
                                                LutSynthParams>::type;

  using SpfdMeth = typename dependency::window::OracleDispatch<
      StaticParams::kSpfdMeth, Ntk, WindowManager, WindowSimulator,
      StaticParams::kMaxCutsSize>::type;

 public:
  AreaResynImpl(Ntk& ntk, Params ps, AreaResynStats& st)
      : ntk_(ntk),
        ps_(ps),
        st_(st),
        wm_params_{ps_.max_num_divisors, ps_.max_reconvergence, ps_.odc_levels},
        window_{WindowManager(ntk_, wm_params_, wm_stats_),
                WindowSimulator(ntk_)},
        synth_(ntk, window_.manager, window_.simulator, st_.st_synth),
        spfds_(ntk, window_.manager, window_.simulator) {}

  void Run();

 private:
  // === Phase decomposition ===
  void SetupTrackers();
  void OptimizeNode(Node n);

  // === Orchestration (delegates to strategy classes) ===
  template <uint32_t NumOutputs>
  bool Orchestrate();

  // === High-level optimization ===
  bool Try0ResubSolo();
  bool TryRewireSolo();
  bool TryStructSolo();
  bool TryWindowSolo();

  // === Shared helpers ===
  std::optional<Signal> Find0ResubSolo();
  template <uint32_t kNumInputs>
  bool TrySynthesis(math::static_vector<Signal, kNumInputs> const&);

 private:
  Ntk& ntk_;
  Params ps_;
  AreaResynStats& st_;
  WindowManagerParams wm_params_;
  windowing::window_manager_stats wm_stats_;
  SharedOptoWindow window_;
  SharedOptoTracks tracks_;
  LutSynth synth_;
  SpfdMeth spfds_;
};

// Implementation of Run split into steps
template <class Ntk, typename Params, typename StaticParams>
void AreaResynImpl<Ntk, Params, StaticParams>::Run() {
  SetupTrackers();

  auto const size_init = ntk_.num_gates();
  ntk_.foreach_gate([&](auto const& n) {
    if (st_.GetNumChanges() >= ps_.stop_after_nth) return false;
    if (n > size_init || ntk_.is_dead(n) || ntk_.is_pi(n)) return true;
    if (ntk_.fanin_size(n) == 1) return true;
    OptimizeNode(n);
    return true;
  });
}

template <class Ntk, typename Params, typename StaticParams>
void AreaResynImpl<Ntk, Params, StaticParams>::SetupTrackers() {
  if (ps_.preserve_times) {
    tracks_.arrivals = new ArrivalsTracker{ntk_};
    double const worst_delay = tracks_.arrivals->worst_delay();
    if (ps_.use_required) {
      tracks_.required = new RequiredTracker{ntk_, worst_delay};
    }
  }
}

template <class Ntk, typename Params, typename StaticParams>
void AreaResynImpl<Ntk, Params, StaticParams>::OptimizeNode(Node n) {
  window_.root = n;

  mockturtle::call_with_stopwatch(st_.time_window_constr,
                                  [&]() { window_.manager.run(n); });

  mockturtle::call_with_stopwatch(st_.time_window_simula, [&]() {
    window_.simulator.run(window_.manager);
  });

  if (ntk_.num_outputs(n) == 1)
    Orchestrate<1>();
  else if (ntk_.num_outputs(n) == 2)
    Orchestrate<2>();
  else {
    // gracefully skip instead of throwing
    fmt::print("[warn] Skipping node {} with {} outputs\n", n,
               ntk_.num_outputs(n));
  }
}

template <class Ntk, typename Params, typename StaticParams>
template <uint32_t NumOutputs>
bool AreaResynImpl<Ntk, Params, StaticParams>::Orchestrate() {
  if constexpr (NumOutputs == 1) {
    if (mockturtle::call_with_stopwatch(st_.time_0resub_solo,
                                        [&]() { return Try0ResubSolo(); }))
      return true;

    if (mockturtle::call_with_stopwatch(st_.time_rewire_solo,
                                        [&]() { return TryRewireSolo(); }))
      return true;

    if (window_.manager.num_mffc() <= 1) return true;

    if (mockturtle::call_with_stopwatch(st_.time_struct_solo,
                                        [&]() { return TryStructSolo(); }))
      return true;

    if (mockturtle::call_with_stopwatch(st_.time_window_solo,
                                        [&]() { return TryWindowSolo(); }))
      return true;

    return true;
  } else if constexpr (NumOutputs == 2) {
    // TODO: dual-output strategy
    return true;
  } else {
    fmt::print("[warn] Unsupported {}-output node\n", NumOutputs);
    return false;
  }
}

template <class Ntk, typename Params, typename StaticParams>
bool AreaResynImpl<Ntk, Params, StaticParams>::Try0ResubSolo() {
  if (!ps_.try_0resub_solo) return false;
  std::optional<Signal> div = Find0ResubSolo();
  if (div) {
    st_.Record0ResubSolo();
    ntk_.substitute_node(window_.root, *div);
    return true;
  }
  return false;
}

template <class Ntk, typename Params, typename StaticParams>
std::optional<typename Ntk::signal>
AreaResynImpl<Ntk, Params, StaticParams>::Find0ResubSolo() {
  Signal froot = ntk_.make_signal(window_.root);
  auto const& signature_root = window_.simulator.Get(froot);
  std::optional<Signal> div;

  window_.manager.foreach_divisor([&](auto const& f, auto i) {
    if (div) return;
    if (ps_.preserve_times) {
      double const time_max = ps_.use_required
                                  ? tracks_.required->get_time(froot)
                                  : tracks_.arrivals->get_time(froot);
      if (tracks_.arrivals->get_time(f) > time_max) return;
    }
    auto const& signature_div = window_.simulator.Get(f);

    auto const diff = signature_div ^ signature_root;
    if (ps_.odc_levels <= 0) {
      auto const num_ones = kitty::count_ones(diff);
      if (num_ones == 0)
        div = f;
      else if (num_ones == diff.num_bits())
        div = ntk_.create_not(f);
    } else {
      auto const& care = window_.simulator.get_careset();
      auto const num_ones = kitty::count_ones(boolean::binary_and(diff, care));
      if (num_ones == 0)
        div = f;
      else if (num_ones == kitty::count_ones(care))
        div = ntk_.create_not(f);
    }
  });
  return div;
}

template <class Ntk, typename Params, typename StaticParams>
bool AreaResynImpl<Ntk, Params, StaticParams>::TryRewireSolo() {
  if (!ps_.try_rewire_solo) return false;

  Signal const froot = ntk_.make_signal(window_.root);
  auto const& care = window_.simulator.get_careset();
  auto const& func = window_.simulator.Get(froot);

  // Base per-fanin simulations by value (no dangling pointers)
  std::vector<Signature> sims_base;
  sims_base.reserve(ntk_.fanin_size(window_.root));
  ntk_.foreach_fanin(window_.root, [&](auto const& fi, auto) {
    sims_base.push_back(window_.simulator.Get(fi));
  });

  std::optional<std::pair<Signal, Signal>> rewire;

  ntk_.foreach_fanin(window_.root, [&](auto const& fi, auto ii) {
    if (rewire) return;

    auto const ni = ntk_.get_node(fi);
    // Candidate removal: non-PI, multi-fanin, ONLY fanout is the root
    bool removable = !ntk_.is_pi(ni) && (ntk_.fanin_size(ni) > 1);
    ntk_.foreach_fanout(
        ni, [&](auto const& no) { removable &= (no == window_.root); });
    if (!removable) return;

    // Create a modified copy with the ii-th input flipped
    auto sims_mod = sims_base;
    ntk_.foreach_fanin(window_.root, [&](auto const& fii, auto iii) {
      if (fii == fi) sims_mod[ii] = ~sims_mod[iii];
    });

    // Build a transient pointer view expected by compute()
    std::vector<Signature const*> sims_view;
    sims_view.reserve(sims_mod.size());
    for (auto& s : sims_mod) sims_view.push_back(&s);

    // Simulate root with flipped input ii
    Signature const sim_flip = ntk_.compute(froot, sims_view);

    // Care for this input: only minterms where flipping ii changes the root
    Signature const care_inp =
        boolean::binary_and(care, boolean::binary_xor(sim_flip, func));

    // Try to match input `fi` with a divisor under care_inp
    auto const& signature_inp = window_.simulator.Get(fi);

    window_.manager.foreach_divisor([&](auto const& d, auto) {
      if (rewire) return;

      if (ps_.preserve_times) {
        double const time_max = ps_.use_required
                                    ? tracks_.required->get_time(fi)
                                    : tracks_.arrivals->get_time(fi);
        if (tracks_.arrivals->get_time(d) > time_max) return;
      }

      auto const& signature_div = window_.simulator.Get(d);
      auto const diff = boolean::binary_xor(signature_div, signature_inp);
      auto const ones = kitty::count_ones(boolean::binary_and(diff, care_inp));

      if (ones == 0) {
        rewire = {fi, d};
      } else if (ones == kitty::count_ones(care_inp)) {
        rewire = {fi, ntk_.create_not(d)};
      }
    });
  });

  if (rewire) {
    // Rebuild root with one child replaced
    auto tt = ntk_.node_function(window_.root);
    auto children = ntk_.get_children(window_.root);
    for (auto& c : children) {
      if (c == rewire->first) {
        c = rewire->second;
      }
    }
    Signal const fnew = ntk_.create_node(children, tt);
    ntk_.substitute_node(window_.root, fnew);
    st_.RecordRewireSolo();
    return true;
  }
  return false;
}

template <class Ntk, typename Params, typename StaticParams>
bool AreaResynImpl<Ntk, Params, StaticParams>::TryStructSolo() {
  if (!ps_.try_struct_solo) return false;

  auto const& inps = window_.manager.get_inputs();

  if (TrySynthesis<WindowManager::inps_capacity>(inps)) {
    st_.RecordStructSolo();
    return true;
  }
  return false;
}

template <class Ntk, typename Params, typename StaticParams>
template <uint32_t kNumInputs>
bool AreaResynImpl<Ntk, Params, StaticParams>::TrySynthesis(
    math::static_vector<Signal, kNumInputs> const& inps) {
  Signal froot = ntk_.make_signal(window_.root);
  Signature const& func = window_.simulator.Get(froot);
  bool const synthesized =
      ps_.preserve_times
          ? synth_.template Run<kNumInputs>(func, inps, *tracks_.arrivals)
          : synth_.template Run<kNumInputs>(func, inps);
  if (!synthesized) return false;
  if (synth_.NumGates() <= window_.manager.num_mffc()) {
    if (ps_.preserve_times) {
      double const time_max = ps_.use_required
                                  ? tracks_.required->get_time(froot)
                                  : tracks_.arrivals->get_time(froot);
      if (synth_.GetTime() > time_max) return false;
    }
    auto fnew = synth_.template Insert<kNumInputs>(inps);
    ntk_.substitute_node(window_.root, fnew);
    return true;
  }
  return false;
}

template <class Ntk, typename Params, typename StaticParams>
bool AreaResynImpl<Ntk, Params, StaticParams>::TryWindowSolo() {
  if (!ps_.try_window_solo) return false;
  spfds_.Init();
  auto cut = spfds_.GetNext();
  if (cut) {
    if (TrySynthesis<StaticParams::kMaxCutsSize>(*cut)) {
      st_.RecordWindowSolo();
      return true;
    }
  }

  return false;
}

}  // namespace detail

// Public API
template <class Ntk, typename DynamicParams = DefaultAreaResynDynamicParams,
          typename StaticParams = DefaultAreaResynStaticParams>
void AreaResyn(Ntk& ntk, DynamicParams ps = {}, AreaResynStats* pst = nullptr) {
  AreaResynStats st;
  detail::AreaResynImpl<Ntk, DynamicParams, StaticParams> p(ntk, ps, st);
  p.Run();
  if (pst) *pst = st;
}

}  // namespace array
}  // namespace algorithms
}  // namespace opto
}  // namespace rinox

#endif  // RINOX_ALGORITHMS_ARRAY_AREA_RESYN_HPP_
