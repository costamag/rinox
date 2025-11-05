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
  \file timing.hpp
  \brief Compute timing information of a network.

  This engine can be used for power analysis of mapped network.
  For each node, the following information is stored:
  - The sensing time : the first time at which a transition can happen
  - The sensing time : the first time at which the output is stable
  - A vector of simulation patterns identifying quantized timesteps in this
  interval The class has the following template parameters: \param N : number of
  timesteps used for simulating each signal \param I : number of inputs to be
  used for simulating the static truth table ( 2^I input pairs )

  \author Andrea Costamagna
*/

#pragma once

#include "../network/signal_map.hpp"

namespace rinox {

namespace windowing {

template <class Ntk, uint32_t CubeSizeLeaves = 12>
class window_simulator {
 public:
  using signal_t = typename Ntk::signal;
  using node_index_t = typename Ntk::node;
  static constexpr uint32_t kCubeSizeLeaves = CubeSizeLeaves;
  using signature_t = kitty::static_truth_table<CubeSizeLeaves>;
  using TruthTable = kitty::static_truth_table<CubeSizeLeaves>;
  static constexpr uint32_t num_bits = 1u << CubeSizeLeaves;

 public:
  window_simulator(Ntk& ntk) : ntk_(ntk), sig_to_sim_(ntk) {
    sims_.reserve(1000u);
    init();
  }

  template <typename WinMngr>
  void run(WinMngr const& window) {
    sims_.reserve(window.size());
    sims_.resize(CubeSizeLeaves);
    sig_to_sim_.reset();
    sig_to_sim_.resize();

    assign_inputs(window);
    
    window.foreach_divisor([&](auto const& f, auto i) {
      auto const n = ntk_.get_node(f);
      if (f.output == 0 && !window.is_input(n)) {
        compute(window, n);
      }
    });

    window.foreach_mffc([&](auto const& n, auto i) { compute(window, n); });

    window.foreach_tfo([&](auto const& n, auto i) { compute(window, n); });

    care_ = compute_observability_careset(window);
  }

  signature_t const& get(signal_t const& f) const {
    return sims_[sig_to_sim_[f]];
  }

  signature_t const& Get(signal_t const& f) const {
    return sims_[sig_to_sim_[f]];
  }

  signature_t const get_careset() const { return care_; }

  template <typename WinMngr>
  signature_t const compute_observability_careset(WinMngr const& window) {
    signature_t care;
    auto n = window.get_pivot();
    auto const& outputs = window.get_outputs();
    if (ntk_.num_outputs(n) == outputs.size()) {
      bool no_odcs = true;
      for (auto i = 0u; i < outputs.size(); ++i) {
        no_odcs &= (ntk_.get_node(outputs[i]) == window.get_pivot());
      }
      if (no_odcs) {
        return ~care;
      }
    }

    std::vector<signature_t> old_sims;
    window.foreach_output([&](auto fo, auto io) {
      (void)io;
      old_sims.push_back(sims_[sig_to_sim_[fo]]);
    });

    for (uint32_t m = 1u; m < (1u << ntk_.num_outputs(n)); ++m) {
      int i = 0;
      ntk_.foreach_output(n, [&](auto const& f) {
        if (((m >> i) & 0x1) > 0)
          sims_[sig_to_sim_[f]] = ~sims_[sig_to_sim_[f]];
        i++;
      });
      window.foreach_tfo([&](auto no, auto io) { re_compute(window, no); });

      window.foreach_output([&](auto fo, auto io) {
        re_compute(window, fo);
        care |= (old_sims[io] ^ sims_[sig_to_sim_[fo]]);
      });

      ntk_.foreach_output(n, [&](auto const& f) {
        sims_[sig_to_sim_[f]] = ~sims_[sig_to_sim_[f]];
      });
      window.foreach_tfo([&](auto no, auto io) { re_compute(window, no); });
    }

    return care;
  }

  template <typename WinMngr>
  void compute(WinMngr const& window, node_index_t const& n) {
    static const signature_t const0 = sims_[0].construct();
    if (window.is_input(n)) return;
    if (ntk_.is_constant(n)) {
      if (ntk_.constant(n)) {
        ntk_.foreach_output(n, [&](auto const& fo) {
          sig_to_sim_[fo] = sims_.size();
          sims_.emplace_back((~const0) ^ const0);
        });
      } else {
        ntk_.foreach_output(n, [&](auto const& fo) {
          sig_to_sim_[fo] = sims_.size();
          sims_.emplace_back(const0 ^ const0);
        });
      }
      return;
    }

    std::vector<signature_t const*> sim_ptrs;
    ntk_.foreach_fanin(n, [&](auto const& fi, auto ii) {
      sim_ptrs.push_back(&sims_[sig_to_sim_[fi]]);
    });

    auto const tts = ntk_.compute(n, sim_ptrs);
    uint32_t io = 0;
    ntk_.foreach_output(n, [&](auto const& fo) {
      sig_to_sim_[fo] = sims_.size();
      sims_.push_back(tts[io++]);
    });
    return;
  }

  template <typename WinMngr>
  void re_compute(WinMngr const& window, node_index_t const& n) {
    if (window.is_input(n)) return;

    std::vector<signature_t const*> sim_ptrs;
    ntk_.foreach_fanin(n, [&](auto const& fi, auto ii) {
      sim_ptrs.push_back(&sims_[sig_to_sim_[fi]]);
    });
    auto const tts = ntk_.compute(n, sim_ptrs);
    uint32_t io = 0;
    ntk_.foreach_output(
        n, [&](auto const& fo) { sims_[sig_to_sim_[fo]] = tts[io++]; 
        });
    return;
  }

  template <typename WinMngr>
  void re_compute(WinMngr const& window, signal_t const& f) {
    if (window.is_input(ntk_.get_node(f))) return;

    std::vector<signature_t const*> sim_ptrs;
    ntk_.foreach_fanin(f, [&](auto const& fi, auto ii) {
      sim_ptrs.push_back(&sims_[sig_to_sim_[fi]]);
    });
    auto const tt = ntk_.compute(f, sim_ptrs);
    sims_[sig_to_sim_[f]] = tt;
    return;
  }

 private:
  void init() {
    sims_.resize(CubeSizeLeaves);
    for (auto i = 0u; i < CubeSizeLeaves; ++i)
      kitty::create_nth_var(sims_[i], i);
  }

  template <typename WinMngr>
  void assign_inputs(WinMngr const& window) {
    window.foreach_input([&](auto const& f, auto i) { sig_to_sim_[f] = i; 
    });
  }

 private:
  Ntk& ntk_;
  std::vector<signature_t> sims_;
  signature_t care_;
  network::incomplete_signal_map<uint32_t, Ntk> sig_to_sim_;
};

}  // namespace windowing

}  // namespace rinox