/* rinox: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
  \file database_mining.hpp
  \brief Mining tool to extract high-quality sub-networks from mapped networks

  \author Andrea Costamagna
*/
#pragma once

#include <filesystem>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/views/cell_view.hpp>

#include "../boolean/boolean.hpp"
#include "../dependency/struct_dependencies.hpp"
#include "../evaluation/chains.hpp"
#include "../network/network.hpp"
#include "../windowing/window_manager_tiered.hpp"
#include "../windowing/window_simulator_tiered.hpp"
#include "mapped_database.hpp"

namespace rinox {

namespace databases {

struct MiningParamsDefault : windowing::WindowManagerParamsDefault {
  static constexpr uint32_t kMaxNumLeaves = 12u;
  static constexpr uint32_t num_vars_sign = 12u;
  static constexpr uint32_t kInpsCapacity = 512u;
  static constexpr uint32_t kDivsCapacity = 100u;
  static constexpr uint32_t kMffcCapacity = 512u;
  static constexpr uint32_t kTfosCapacity = 512u;
  static constexpr uint32_t kOutsCapacity = 512u;
  static constexpr uint32_t kMaxCutsSize = 6u;
  static constexpr uint32_t max_cuts_size = 6u;
  static constexpr uint32_t kMaxCutReconv = 3u;
  static constexpr uint32_t max_reconv = 3u;
  static constexpr uint32_t kMaxNumVarsDb = 6u;
  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 0;
  int max_num_tiers_tfi = 10;
  int max_num_tiers_tfo = 0;
  bool delay_constrained = false;
};

struct MiningStats : windowing::WindowManagerStatsDefault {
  uint32_t num_added = 0;

  void Report() {
    std::cout << "Number of added sub-networks: " << num_added << std::endl;
  }
};

template <typename Ntk, typename Params = MiningParamsDefault>
class DatabaseMining {
  static constexpr network::design_type_t Design = Ntk::design_t;
  using Database = mapped_database<Ntk, Params::kMaxNumVarsDb>;
  using Function = kitty::static_truth_table<Params::kMaxNumVarsDb>;
  using Gate = mockturtle::gate;
  using Library = libraries::augmented_library<Design>;
  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using StructDependencies = dependency::struct_dependencies<Ntk, Params>;
  using WindowManagerT = windowing::WindowManager<Ntk, Params>;
  using WindowSimulatorT =
      windowing::WindowSimulator<Ntk, Params::kMaxNumLeaves>;
  using Chain = typename evaluation::chains::bound_chain<Design>;
  using Stats = MiningStats;

 public:
  DatabaseMining(std::vector<Gate> const&, Params const&, Stats&);
  void Collect(Ntk&);

 private:
  std::vector<Gate> gates_;
  Library library_;
  Database db_;
  Params const& ps_;
  Stats& st_;
};

template <typename Ntk, typename Params>
DatabaseMining<Ntk, Params>::DatabaseMining(std::vector<Gate> const& gates,
                                            Params const& ps, Stats& st)
    : library_(gates), db_(library_), ps_(ps), st_(st) {}

template <typename Ntk, typename Params>
void DatabaseMining<Ntk, Params>::Collect(Ntk& ntk) {
  StructDependencies dependencies(ntk);

  WindowManagerT window_manager(ntk, ps_, st_);
  WindowSimulatorT window_simulator(ntk);

  ntk.foreach_gate([&](Node const& n) {
    window_manager.Run(n);
    if (!window_manager.IsValid()) return;
    window_simulator.Run(window_manager);
    dependencies.run(window_manager, window_simulator);
    dependencies.foreach_cut([&](auto& cut, auto) {
      std::vector<Signal> inputs;
      for (auto const& leaf : cut.leaves) {
        inputs.push_back(leaf);
      }
      if (db_.add(ntk, inputs, ntk.make_signal(n))) {
        st_.num_added++;
      }
    });
  });
  st_.Report();
}

}  // namespace databases

}  // namespace rinox