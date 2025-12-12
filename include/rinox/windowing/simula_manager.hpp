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
  \file window_manager.hpp
  \brief Construction of windows for mapped networks.

  \author Andrea Costamagna
*/

#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <mockturtle/utils/node_map.hpp>
#include <optional>
#include <queue>
#include <stack>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../math/math.hpp"
#include "../network/visualize.hpp"

namespace rinox {

namespace windowing {

struct DefaultSimulaManagerParams {
  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 10;
  int32_t num_levels_tfo = 0;
  static constexpr uint32_t kMaxNumLeaves = 20u;
  static constexpr uint32_t kNumVarsSigns = 12u;
};

struct SimulaManagerStats {
  bool valid = false;
};

template <class Ntk, typename Params = DefaultSimulaManagerParams>
class SimulaManager {
 public:
  struct WindowManagerParams {
    uint32_t max_num_divisors;
    uint32_t max_reconvergence;
    int32_t num_levels_tfo;
    static constexpr uint32_t max_num_leaves = Params::kMaxNumLeaves;
  };

  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using WindowManager = window_manager<Ntk, WindowManagerParams>;
  using Stats = SimulaManagerStats;
  static constexpr uint32_t kNumVarsSigns = Params::kNumVarsSigns;
  static constexpr uint32_t kMaxNumLeaves = Params::kMaxNumLeaves;
  using Signature = kitty::static_truth_table<kNumVarsSigns>;

  SimulaManager(Ntk& ntk, Params& ps, Stats& st)
      : ntk_(ntk), sig_to_sim_(ntk), ps_(ps), st_(st){};

 private:
  Ntk& ntk_;
  Params const& ps_;
  Stats& st_;
  std::vector<Signature> sims_;
  network::incomplete_signal_map<uint32_t, Ntk> sig_to_sim_;
};

}  // namespace windowing
}  // namespace rinox
