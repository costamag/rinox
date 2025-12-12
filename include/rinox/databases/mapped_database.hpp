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
  \file mapped_database.hpp
  \brief Manager for databases of mapped networks

  \author Andrea Costamagna
*/
#pragma once

#include <kitty/kitty.hpp>

#include "../boolean/boolean.hpp"
#include "../evaluation/evaluation.hpp"
#include "../io/utils/reader.hpp"
#include "../io/verilog/verilog.hpp"
#include "../io/verilog/write_verilog.hpp"
#include "../math/math.hpp"
#include "../network/utils.hpp"

namespace rinox {

namespace databases {

/*! \brief Database of mapped networks
 *
 * \tparam NtkDb Network type of the stored database
 */
template <typename NtkDb, uint32_t MaxNumVars = 6u>
class mapped_database {
 public:
  static constexpr network::design_type_t design_t = NtkDb::design_t;
  static constexpr uint32_t max_num_vars = MaxNumVars;
  using chain_t = typename evaluation::chains::bound_chain<design_t>;
  using signal_t = typename NtkDb::signal;
  using signals_t = rinox::math::static_vector<signal_t, MaxNumVars>;
  using library_t = typename libraries::augmented_library<design_t>;
  using node_index_t = typename NtkDb::node_index_t;
  using truth_table_t = kitty::static_truth_table<MaxNumVars>;
  using chain_simulator_t = evaluation::chain_simulator<chain_t, truth_table_t>;

  struct database_entry_t {
    bool operator<(database_entry_t const& other) const {
      bool dominates = true;
      dominates &= area < other.area;
      dominates &= switches < other.switches;
      bool one_strict = false;
      for (int i = 0; i < delays.size(); ++i) {
        one_strict |= delays[i] < other.delays[i];
        dominates &= delays[i] <= other.delays[i];
      }
      return dominates & one_strict;
    }

    bool operator>=(database_entry_t const& other) const {
      bool dominated = true;
      dominated &= area >= other.area;
      dominated &= switches >= other.switches;
      for (int i = 0; i < delays.size(); ++i) {
        dominated &= delays[i] >= other.delays[i];
      }
      return dominated;
    }

    /*! \brief Area of the sub-network */
    double area;
    /*! \brief Zero-delay switching */
    uint32_t switches;
    /*! \brief Longest path from each pin to the outputs */
    std::vector<double> delays;
    /*! \brief Node implementing the functionality */
    node_index_t index;
  };

  struct database_row_t {
    size_t size() const { return entries.size(); }

    void push_back(database_entry_t const& entry) { entries.push_back(entry); }

    database_entry_t& operator[](uint32_t i) { return entries[i]; }

    boolean::symmetries_t symm;
    truth_table_t repr;
    std::vector<database_entry_t> entries;
  };

  struct match_t {
    match_t() = default;
    match_t(match_t&&) = default;
    match_t(match_t const&) = default;

    match_t& operator=(match_t const&) = default;
    match_t& operator=(match_t&&) = default;

    match_t(boolean::permutation_t const& perm, uint64_t const& row)
        : perm(perm), row(row) {}

    boolean::permutation_t perm;
    uint64_t row;
  };

 public:
  mapped_database(library_t& lib) : lib_(lib), ntk_(lib), simulator_(lib) {
    for (auto i = 0u; i < MaxNumVars; ++i) pis_.push_back(ntk_.create_pi());

    for (int i = 0; i < MaxNumVars; ++i) {
      kitty::create_nth_var(proj_funcs_[i], i);
    }
    /* extract the chain's functionality */
    for (auto i = 0u; i < MaxNumVars; ++i) {
      sims_ptrs_.push_back(&proj_funcs_[i]);
    }
    AddGates();
  }

  mapped_database() = default;
  mapped_database(const mapped_database&) = default;
  mapped_database(mapped_database&&) = default;
  mapped_database& operator=(const mapped_database&) = default;
  mapped_database& operator=(mapped_database&&) = default;
  ~mapped_database() = default;

  library_t& get_library() { return lib_; }

#pragma region Init
  void AddGates() {
    NtkDb ntk_tmp(lib_);
    rinox::math::static_vector<signal_t, MaxNumVars> pis;

    for (auto i = 0u; i < MaxNumVars; ++i) pis.push_back(ntk_tmp.create_pi());

    auto const& raw_gates = lib_.get_raw_gates();
    for (auto const& g : raw_gates) {
      if (g.num_vars <= MaxNumVars) {
        std::vector<signal_t> children;
        for (int i = 0; i < g.num_vars; ++i) children.push_back(pis[i]);
        auto f = ntk_tmp.create_node(children, g.id);
        add(ntk_tmp, children, f);
      }
    }
  }
#pragma endregion

#pragma region loading

  void load(std::string const& filename) {
    std::ifstream in(lorina::detail::word_exp_filename(filename),
                     std::ifstream::in);
    load(in);
  }

  void load(std::istream& file) {
    /* load the database */
    auto const& gates = lib_.get_raw_gates();
    NtkDb ntk_ext(gates);
    rinox::io::reader<NtkDb> read(ntk_ext);
    auto parse_status = rinox::io::verilog::read_verilog(file, read);
    assert(parse_status == lorina::return_code::success);
    /* output-by-output add the entry */
    std::vector<typename NtkDb::signal> inputs;
    ntk_ext.foreach_pi(
        [&](auto const& n) { inputs.push_back(ntk_ext.make_signal(n)); });
    ntk_ext.foreach_po([&](auto const& f) { add(ntk_ext, inputs, f); });
  }

#pragma endregion

#pragma region Saving

  void commit(std::string const& file) {
    io::verilog::write_verilog(ntk_, file);
  }

  void commit(std::ostream& os) { io::verilog::write_verilog(ntk_, os); }

#pragma endregion

#pragma region Getters
  /*! \brief Get the number of rows in the database */
  uint64_t num_rows() const { return database_.size(); }

  /*! \brief Get the number of sub-networks stored */
  uint64_t size() const { return ntk_.num_pos(); }
#pragma endregion

#pragma region Insert in Database
  uint64_t memoize_func(truth_table_t const& tt) {
    uint64_t row_index;
    const auto it0 = func_to_match_.find(tt);
    if (it0 != func_to_match_.end()) {
      row_index = (it0->second).row;
    } else {
      auto [repr, _, perm] = kitty::exact_p_canonization(tt);
      const auto it1 = repr_to_row_.find(repr);
      /* insert in the database */
      if (it1 != repr_to_row_.end()) {
        row_index = it1->second;
      } else  // the representative doesn't have any implementation yet
      {
        row_index = database_.size();
        database_.emplace_back();
        database_.back().symm = boolean::symmetries_t(repr);
        database_.back().repr = repr;
        repr_to_row_[repr] = row_index;
      }

      func_to_match_[tt] = match_t{boolean::permutation_t(perm), row_index};
    }
    return row_index;
  }

  /*! \brief Insert a mapped chain into the database
   *
   * This method inserts the chain into the database the way it is provided.
   * Any check for whether the chain should be inserted or not should be handled
   * at a higher level of abstraction.
   */
  bool add(evaluation::chains::bound_chain<design_t> chain) {
    simulator_(chain, sims_ptrs_);
    truth_table_t const tt =
        simulator_.get_simulation(chain, sims_ptrs_, chain.po_at(0));
    /* perform P-canonization on the chain */
    uint64_t row = memoize_func(tt);

    rinox::evaluation::chains::perm_canonize(chain, func_to_match_[tt].perm);
    rinox::evaluation::chains::time_canonize(chain, lib_, database_[row].symm);

    bool is_inserted = add(chain, row);

    simulator_(chain, sims_ptrs_);
    truth_table_t const tt2 =
        simulator_.get_simulation(chain, sims_ptrs_, chain.po_at(0));
    assert(kitty::equal(tt2, database_[row].repr));
    /* create the new entry */

    return is_inserted;
  }

  /*! \brief Insert a mapped sub-network into the database
   *
   * The sub-network is identified by the input and output signals. The
   * subnetwork is first converted to a chain, which is then inserted into the
   * database.
   */
  template <typename Ntk>
  bool add(Ntk& ntk, std::vector<typename Ntk::signal> const& inputs,
           typename Ntk::signal const& output) {
    evaluation::chains::bound_chain<design_t> chain(MaxNumVars);
    rinox::evaluation::chains::extract(chain, ntk, inputs, output);
    return add(chain);
  }

 private:
  bool add(evaluation::chains::bound_chain<design_t>& chain, uint64_t row) {
    database_entry_t entry;
    entry.area = chain.get_area(lib_);
    entry.switches = simulator_.get_switches(chain);
    entry.delays = get_longest_paths(chain, lib_);

    for (auto i = 0u; i < database_[row].size(); ++i) {
      if (entry >= database_[row][i]) {
        return false;
      } else if (entry < database_[row][i]) {
        auto const f =
            rinox::evaluation::chains::insert<NtkDb, NtkDb::design_t, signals_t,
                                              true>(ntk_, pis_, chain);
        ntk_.substitute_node(database_[row][i].index, f);
        entry.index = ntk_.get_node(f);
        database_[row][i] = entry;
        return true;
      }
    }
    // TODO: add capacity
    auto const f =
        rinox::evaluation::chains::insert<NtkDb, NtkDb::design_t, signals_t,
                                          true>(ntk_, pis_, chain);
    if (ntk_.is_po(f)) return false;  // do not re-insert POs in the database
    ntk_.create_po(f);
    entry.index = ntk_.get_node(f);
    database_[row].push_back(entry);
    return true;
  }
#pragma endregion

#pragma region Lookup
 public:
  template <typename Times, typename... Vecs>
  std::optional<uint64_t> boolean_matching(truth_table_t const& func,
                                           Times& times, Vecs&... others) {
    // Resize times
    times.resize(MaxNumVars, std::numeric_limits<double>::max());

    // Resize all other vectors
    (((void)others.resize(MaxNumVars)), ...);

    // check if the representative is already in the database
    auto match = get_match(func);
    if (match) {
      auto const row_index = (*match).row;
      database_row_t const& row = database_[row_index];
      auto const perm = (*match).perm;
      auto const symm = row.symm;

      // Apply permutation to all vectors simultaneously
      boolean::forward_permute_inplace(perm, times, others...);

      // Apply time-based symmetric sorting to all vectors simultaneously
      boolean::sort_symmetric(
          symm, [&](auto const& a, auto const& b) { return a < b; }, times,
          others...);

      return row_index;
    }
    return std::nullopt;
  }

  template <typename Fn>
  void foreach_entry(uint64_t row_index, Fn&& fn) {
    for (database_entry_t const& entry : database_[row_index].entries) {
      fn(entry);
    }
  }

  template <typename Ntk, typename Signals>
  node_index_t write(database_entry_t const& entry, Ntk& ntk,
                     Signals const& leaves) {
    return write(entry.index, ntk, leaves);
  }

  template <typename Ntk, typename Signals>
  node_index_t write(typename NtkDb::node const& index, Ntk& ntk,
                     Signals const& leaves) {
    ntk_.incr_trav_id();

    std::function<node_index_t(typename NtkDb::node const&)> insert;

    insert = [&](typename NtkDb::node const& n) -> node_index_t {
      if (ntk_.visited(n) == ntk_.trav_id()) return ntk_.value(n);

      if (ntk_.is_pi(n)) {
        ntk_.set_value(n, ntk.get_node(leaves[ntk_.pi_index(n)]));
        ntk_.set_visited(n, ntk_.trav_id());
        return ntk_.value(n);
      }

      std::vector<signal_t> children(ntk_.fanin_size(n));
      ntk_.foreach_fanin(n, [&](auto const& fi, auto i) {
        auto const ni = ntk_.get_node(fi);
        auto const nnew = insert(ni);
        children[i] = signal_t{nnew, fi.output};
      });

      auto const ids = ntk_.get_binding_ids(n);
      auto nnew = ntk.get_node(ntk.create_node(children, ids));
      ntk_.set_value(n, nnew);
      ntk_.set_visited(n, ntk_.trav_id());

      return nnew;
    };

    return insert(index);
  }

 private:
  std::optional<match_t> get_match(truth_table_t const& tt) {
    const auto it0 = func_to_match_.find(tt);
    if (it0 != func_to_match_.end()) {
      return (it0->second);
    } else {
      auto [repr, _, perm] = kitty::exact_p_canonization(tt);
      const auto it1 = repr_to_row_.find(repr);
      /* insert in the database */
      if (it1 != repr_to_row_.end()) {
        uint64_t const row = it1->second;
        auto const match = match_t{boolean::permutation_t(perm), row};
        func_to_match_[tt] = match;
        return match;
      }
    }
    return std::nullopt;
  }

  template <typename E, typename T>
  void perm_matching(std::vector<E>& leaves, std::vector<T>& times,
                     boolean::permutation_t const& perm) {
    forward_permute_inplace(perm, leaves, times);
  }

  /*! \brief Permutes the input variables to have the ones closest to the output
   * last
   *
   *    3
   *   2
   *  1
   * 0
   * */
  template <typename E, typename T>
  void time_matching(std::vector<E>& leaves, std::vector<T>& times,
                     boolean::symmetries_t const& symm) {
    // sort the symmetric input variables
    sort_symmetric(
        times, symm, [&](auto const& a, auto const& b) { return a < b; },
        leaves);
  }
#pragma endregion

 private:
  /*! \brief Map a truth table to a storage of nodes and input permutations */
  std::vector<database_row_t> database_;

  /*! \brief Map a completely specified truth table to a database row */
  phmap::flat_hash_map<truth_table_t, match_t, kitty::hash<truth_table_t>>
      func_to_match_;
  phmap::flat_hash_map<truth_table_t, uint64_t, kitty::hash<truth_table_t>>
      repr_to_row_;

  /*! \brief Database represented as a network */
  NtkDb ntk_;
  rinox::math::static_vector<signal_t, MaxNumVars> pis_;

  /*! \brief Technology library */
  library_t lib_;

  /*! \brief Simulation engine for mapped chains */
  chain_simulator_t simulator_;

  /*! \brief Vector of simulationpatterns */
  std::array<truth_table_t, MaxNumVars> proj_funcs_;
  std::vector<truth_table_t const*> sims_ptrs_;
};

} /* namespace databases */

} /* namespace rinox */