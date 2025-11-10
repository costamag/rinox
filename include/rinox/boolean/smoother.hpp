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
  \file bdd_synthesizer.hpp
  \brief Synthesis engine from Binary Decision Trees

  \author Andrea Costamagna
*/

#ifndef RINOX_SMOOTHER_HPP_
#define RINOX_SMOOTHER_HPP_

#include "../math/bdd.hpp"
#include "metrics.hpp"

namespace rinox {

namespace boolean {

template <typename CSFunction, typename Var, uint32_t kMaxNumVars>
class Smoother {
 public:
  using ISFunction = kitty::ternary_truth_table<CSFunction>;
  using Simulations = math::static_vector<CSFunction, (1u << kMaxNumVars) + 2u>;
  using Vars = math::static_vector<Var, kMaxNumVars>;
  using SimPtrs = math::static_vector<CSFunction const*, kMaxNumVars>;
  using BDD = math::BinaryDecisionTree<Var, 1u << kMaxNumVars>;

  Smoother() {
    sims_.emplace_back();
    sims_.push_back(~sims_[0]);
  }

  bool IsConst0(CSFunction const& bits, CSFunction const& care) const {
    return kitty::count_ones(bits & care) == 0;
  }

  bool IsConst1(CSFunction const& bits, CSFunction const& care) const {
    return kitty::equal(bits & care, care);
  }

  uint32_t RunCARTRecursive(Vars const& vars, ISFunction const& func) {
    return RunCARTRecursive(vars, func._bits, func._care);
  }

  uint32_t RunCARTRecursive(Vars const& vars, CSFunction const& bits,
                            CSFunction const& care) {
    if (IsConst0(bits, care))
      return BDD::kConst0;
    else if (IsConst1(bits, care))
      return BDD::kConst1;

    std::optional<Var> best;

    int best_info = boolean::GiniImpurity(bits, care);
    for (auto const v : vars) {
      auto const& sim_var = *cache_[v];
      auto const care_else = care & ~sim_var;
      auto const care_then = care & sim_var;
      uint32_t const info = boolean::GiniImpurity(bits, care_else) +
                            boolean::GiniImpurity(bits, care_then);

      if (info < best_info) {
        best_info = info;
        best = v;
      }
    }
    if (!best) throw std::runtime_error("no best variable found");

    auto const& sim_var = *cache_[*best];
    auto const& care_then = care & sim_var;
    auto const& care_else = care & ~sim_var;
    uint32_t t;
    uint32_t e;
    if (IsConst0(bits, care_then))
      t = BDD::kConst0;
    else if (IsConst1(bits, care_then))
      t = BDD::kConst1;
    else
      t = RunCARTRecursive(vars, bits, care_then);

    if (IsConst0(bits, care_else))
      e = BDD::kConst0;
    else if (IsConst1(bits, care_else))
      e = BDD::kConst1;
    else
      e = RunCARTRecursive(vars, bits, care_else);

    return bdd_.AddNode(*best, t, e);
  }

  uint32_t RunFastRecursive(Vars const& vars, ISFunction const& func,
                            uint32_t var_index) {
    return RunFastRecursive(vars, func._bits, func._care, var_index);
  }

  uint32_t RunFastRecursive(Vars const& vars, CSFunction const& bits,
                            CSFunction const& care, uint32_t var_index) {
    if (IsConst0(bits, care))
      return BDD::kConst0;
    else if (IsConst1(bits, care))
      return BDD::kConst1;

    auto const& sim_var = *cache_[vars[var_index]];
    auto const& care_then = care & sim_var;
    auto const& care_else = care & ~sim_var;
    uint32_t t;
    uint32_t e;
    if (IsConst0(bits, care_then))
      t = BDD::kConst0;
    else if (IsConst1(bits, care_then))
      t = BDD::kConst1;
    else
      t = RunFastRecursive(vars, bits, care_then, var_index + 1);

    if (IsConst0(bits, care_else))
      e = BDD::kConst0;
    else if (IsConst1(bits, care_else))
      e = BDD::kConst1;
    else
      e = RunFastRecursive(vars, bits, care_else, var_index + 1);

    return bdd_.AddNode(vars[var_index], t, e);
  }

  CSFunction RunCART(CSFunction const& bits, CSFunction const& care,
                     Vars const& vars, SimPtrs const& sim_ptrs) {
    if (vars.size() != sim_ptrs.size())
      throw std::invalid_argument("variables and sim_ptrs size mismatch");

    cache_.clear();
    for (uint32_t i = 0; i < vars.size(); ++i) cache_[vars[i]] = sim_ptrs[i];

    bdd_.Reset();
    RunCARTRecursive(vars, bits, care);

    return Simulate();
  }

  CSFunction RunCART(ISFunction const& func, Vars const& vars,
                     SimPtrs const& sim_ptrs) {
    return RunCART(func._bits, func._care, vars, sim_ptrs);
  }

  uint32_t RunC4d5Recursive(Vars const& vars, CSFunction const& bits,
                            CSFunction const& care) {
    if (IsConst0(bits, care)) return BDD::kConst0;
    if (IsConst1(bits, care)) return BDD::kConst1;

    std::optional<Var> best;
    double best_info = -1.0;

    for (auto const v : vars) {
      auto const& sim = *cache_.at(v);  // use at() to catch missing
      const uint32_t num_cares_then = kitty::count_ones(care & sim);
      const uint32_t num_cares_else = kitty::count_ones(care & ~sim);

      // Skip non-partitioning variables
      if (num_cares_then == 0 || num_cares_else == 0) continue;

      double const info = boolean::MutualInformation(sim, bits, care);
      if (info > best_info) {
        best_info = info;
        best = v;
      }
    }

    if (!best) {
      // No splitter can reduce care: make a majority leaf (prevents infinite
      // recursion)
      const uint32_t ones = kitty::count_ones(bits & care);
      const uint32_t N = kitty::count_ones(care);
      return (ones * 2 >= N) ? BDD::kConst1 : BDD::kConst0;
    }

    auto const& sim_var = *cache_.at(*best);
    auto const care_then = care & sim_var;
    auto const care_else = care & ~sim_var;

    // Optional but safer: remove 'best' from the var list for children
    Vars child_vars = vars;
    child_vars.erase(std::remove(child_vars.begin(), child_vars.end(), *best),
                     child_vars.end());

    uint32_t t = IsConst0(bits, care_then) ? BDD::kConst0
                 : IsConst1(bits, care_then)
                     ? BDD::kConst1
                     : RunC4d5Recursive(child_vars, bits, care_then);

    uint32_t e = IsConst0(bits, care_else) ? BDD::kConst0
                 : IsConst1(bits, care_else)
                     ? BDD::kConst1
                     : RunC4d5Recursive(child_vars, bits, care_else);

    return bdd_.AddNode(*best, t, e);
  }

  CSFunction RunC4d5(CSFunction const& bits, CSFunction const& care,
                     Vars const& vars, SimPtrs const& sim_ptrs) {
    if (vars.size() != sim_ptrs.size())
      throw std::invalid_argument("variables and sim_ptrs size mismatch");

    cache_.clear();
    for (uint32_t i = 0; i < vars.size(); ++i) cache_[vars[i]] = sim_ptrs[i];

    bdd_.Reset();
    RunC4d5Recursive(vars, bits, care);

    return Simulate();
  }

  CSFunction RunC4d5(ISFunction const& func, Vars const& vars,
                     SimPtrs const& sim_ptrs) {
    return RunC4d5(func._bits, func._care, vars, sim_ptrs);
  }

  CSFunction RunFast(ISFunction const& func, Vars const& vars,
                     SimPtrs const& sim_ptrs) {
    if (vars.size() != sim_ptrs.size())
      throw std::invalid_argument("variables and sim_ptrs size mismatch");

    cache_.clear();
    for (uint32_t i = 0; i < vars.size(); ++i) cache_[vars[i]] = sim_ptrs[i];

    bdd_.Reset();
    RunFastRecursive(vars, func, 0u);

    return Simulate();
  }

  CSFunction Simulate() {
    sims_.resize(2u);
    bdd_.ForeachNode([&](auto const& n, auto const& index) {
      auto const& sim_var = *cache_[n.v];
      auto const& sim_then = sims_[n.t];
      auto const& sim_else = sims_[n.e];
      sims_.push_back((sim_var & sim_then) | (~sim_var & sim_else));
    });
    return sims_[sims_.size() - 1];
  }

  CSFunction RunLabelPropagation(ISFunction const& func) {
    return RunLabelPropagation(func._bits, func._care);
  }

  CSFunction RunLabelPropagation(CSFunction const& bits,
                                 CSFunction const& care) {
    const uint32_t n = care.num_vars();

    std::vector<CSFunction> flips;
    std::vector<uint32_t> dont_cares;
    flips.reserve(care.num_bits());
    dont_cares.reserve(care.num_bits());

    // Precompute neighbor masks for each don't-care minterm
    for (uint32_t m = 0; m < care.num_bits(); ++m) {
      if (!kitty::get_bit(care, m)) {
        dont_cares.push_back(m);
        CSFunction flip{};  // zero-initialize!
        for (uint32_t v = 0; v < n; ++v) {
          kitty::set_bit(flip, m ^ (1u << v));  // Hamming-1 neighbors
        }
        flips.push_back(flip);
      }
    }

    // Initialize: keep cares fixed, randomize don't-cares
    CSFunction new_bits;
    kitty::create_random(new_bits);
    new_bits = (bits & care) | (new_bits & ~care);

    bool changed;
    do {
      changed = false;

      // Asynchronous (in-place) updates help avoid 2-cycles on bipartite graphs
      for (size_t i = 0; i < dont_cares.size(); ++i) {
        const uint32_t m = dont_cares[i];
        const auto& nb = flips[i];

        const uint64_t ones_nb = kitty::count_ones(nb & new_bits);
        const bool next = (ones_nb > n / 2);  // strict majority
        const bool current = kitty::get_bit(new_bits, m);

        // Tie-break: keep current (swap to `false` to bias 0 on ties)
        const bool next_with_tie = (ones_nb * 2 == n) ? current : next;

        if (current != next_with_tie) {
          if (next_with_tie)
            kitty::set_bit(new_bits, m);
          else
            kitty::clear_bit(new_bits, m);
          changed = true;
        }
      }

    } while (changed);

    return new_bits;
  }

 private:
  BDD bdd_;
  Simulations sims_;
  std::unordered_map<Var, CSFunction const*> cache_;
};

}  // namespace boolean

}  // namespace rinox

#endif