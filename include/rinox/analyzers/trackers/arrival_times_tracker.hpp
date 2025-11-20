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
 * \file arrival_times_tracker.hpp
 * \brief Compute arrival time information of a network and updates upon change
 *
 * \author Andrea Costamagna
 */

#pragma once

#include <limits>
#include <mockturtle/networks/events.hpp>

#include "../../network/signal_map.hpp"
#include "../../network/tfo_manager.hpp"

namespace rinox {

namespace analyzers {

namespace trackers {

/*! \brief Engine to evaluate the arrival times of a network.
 *
 * This engine computes the arrival times of a network and keeps them
 up-to-date.
 * During construction it is possible to specify a vector of arrival times. If
 * not provided, the engine assumes 0 arrival time at all the PIs. At
 construction,
 * the arrival times are propagated in the network.
 *
 * The engine is equipped with a transitive fanout (TFO) manager to handle the
 * timing updates. Two events trigger the update:
 * - Node addition: The arrival time of the node is computed from the fanins
 * - Node modification: The arrival time of the TFO of the affected nodes is
 updated.
 *
 * \tparam Ntk the network type to be analyzed.
 *
 * \verbatim embed:rst

   Example

   .. code-block:: c++

      bound_network ntk( gates );
      auto const a = ntk.create_pi();
      auto const f1 = ntk.create_node( { a }, 0 );
      arrival_times_tracker arrival( ntk );
      auto const f2 = ntk.create_node( { f1 }, 0 );
   \endverbatim
 */
template <class Ntk>
class arrival_times_tracker {
 public:
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;

 public:
  arrival_times_tracker(Ntk& ntk) : ntk_(ntk), times_(ntk), tfo_(ntk) {
    init();
  }

  arrival_times_tracker(Ntk& ntk, std::vector<double> const& input_arrivals)
      : ntk_(ntk), times_(ntk), input_(input_arrivals), tfo_(ntk) {
    init();
  }

  void init() {
    /* check if the network implements the needed functionalities */
    static_assert(mockturtle::is_network_type_v<Ntk>,
                  "Ntk is not a network type");
    static_assert(mockturtle::has_incr_trav_id_v<Ntk>,
                  "Ntk does not implement the incr_trav_id method");
    static_assert(mockturtle::has_foreach_fanin_v<Ntk>,
                  "Ntk does not implement the foreach_fanin method");
    static_assert(rinox::traits::has_foreach_output_v<Ntk>,
                  "Ntk does not implement the foreach_output method");
    static_assert(mockturtle::has_foreach_fanout_v<Ntk>,
                  "Ntk does not implement the foreach_fanin method");
    static_assert(mockturtle::has_foreach_pi_v<Ntk>,
                  "Ntk does not implement the foreach_pi method");
    static_assert(mockturtle::has_foreach_po_v<Ntk>,
                  "Ntk does not implement the foreach_pi method");
    static_assert(mockturtle::has_set_value_v<Ntk>,
                  "Ntk does not implement the set_visited method");
    static_assert(mockturtle::has_value_v<Ntk>,
                  "Ntk does not implement the set_visited method");
    static_assert(mockturtle::has_set_visited_v<Ntk>,
                  "Ntk does not implement the set_visited method");
    static_assert(mockturtle::has_visited_v<Ntk>,
                  "Ntk does not implement the set_visited method");
    static_assert(mockturtle::has_get_node_v<Ntk>,
                  "Ntk does not implement the get_node method");
    static_assert(mockturtle::has_size_v<Ntk>,
                  "Ntk does not implement the size method");
    compute_arrival_times();

    add_event_ = ntk_.events().register_add_event([&](const node_index_t& n) {
      times_.resize();
      compute_arrival_time(n);
    });

    modified_event_ = ntk_.events().register_modified_event(
        [&](const auto& n, auto old_children) {
          times_.resize();
          compute_arrival_times();
          //      update_arrival_times_tfo( n );
          //
          //      for ( auto const& f : old_children )
          //      {
          //        update_arrival_times_tfo( ntk_.get_node( f ) );
          //      }
        });
  }

  ~arrival_times_tracker() {
    if (add_event_) {
      ntk_.events().release_add_event(add_event_);
    }

    if (modified_event_) {
      ntk_.events().release_modified_event(modified_event_);
    }
  }

#pragma region Interface methods
 public:
  [[nodiscard]] double worst_delay() const {
    double worst_delay = 0;
    ntk_.foreach_po([&](signal_t const& f) {
      worst_delay = std::max(worst_delay, times_[f]);
    });
    return worst_delay;
  }

  [[nodiscard]] double get_time(signal_t const f) const { return times_[f]; }

  [[nodiscard]] double GetArrival(signal_t const f) const { return times_[f]; }

#pragma endregion

#pragma region Implementation details
 private:
  bool is_marked_ready(node_index_t const& n) {
    return ntk_.value(n) == ntk_.trav_id();
  };

  void make_ready(node_index_t const& n) { ntk_.set_value(n, ntk_.trav_id()); };

  void compute_arrival_times() {
    times_.reset();
    if (input_.size() < ntk_.num_pis()) {
      input_ = std::vector<double>(ntk_.num_pis(), 0);
    }

    ntk_.incr_trav_id();
    ntk_.foreach_pi([&](auto n, auto index) {
      times_[ntk_.make_signal(n)] = input_[index];
      make_ready(n);
    });

    ntk_.foreach_po([&](auto f) { compute_arrival_times_tfi(f); });
  }

  /*! \brief Compute the arrival times of the nodes in the TFI of a signal's
   * node*/
  void compute_arrival_times_tfi(signal_t const& f) {
    node_index_t n = ntk_.get_node(f);
    if (is_marked_ready(n) || ntk_.is_pi(n)) {
      return;
    }
    make_ready(n);

    ntk_.foreach_fanin(
        n, [&](auto fi, auto ii) { compute_arrival_times_tfi(fi); });
    compute_arrival_time(n);
  }

  /*! \brief Efficient update of the arrival times in the TFO of a node.\
   */
  void update_arrival_times_tfo(node_index_t const& n) {
    /* mark all the nodes that must be considered */
    tfo_.init(n);

    std::vector<node_index_t> old_worklist{n};
    std::vector<node_index_t> new_worklist;
    old_worklist.reserve(100);
    new_worklist.reserve(100);
    tfo_.mark_seen(n);

    bool progress{true};
    /* iterate over all the nodes having some inputs ready */
    while (progress && !old_worklist.empty()) {
      progress = false;
      for (node_index_t const& u : old_worklist) {
        bool ready = true;
        ntk_.foreach_fanin(u, [&](auto f, auto i) {
          node_index_t ni = ntk_.get_node(f);
          bool const cond1 = ntk_.is_pi(ni);
          bool const cond2 = !tfo_.belongs_to_tfo(ni);
          bool const cond3 =
              (tfo_.belongs_to_tfo(ni) && tfo_.is_marked_ready(ni));
          ready &= (cond1 || cond2 || cond3);
        });
        tfo_.mark_seen(u);
        if (ready) {
          progress = true;
          tfo_.mark_ready(u);
          ntk_.foreach_output(u, [&](auto const& fu) {
            double old_arrival = times_[fu];
            compute_arrival_time_at_pin(fu);
            ntk_.foreach_fanout(fu, [&](node_index_t const& o) {
              if ((tfo_.belongs_to_tfo(o) &&
                   !tfo_.is_marked_seen(o)))  // only insert if not seen
              {
                new_worklist.push_back(o);
                tfo_.mark_seen(o);
              }
            });
          });
        } else {
          new_worklist.push_back(u);
        }
      }
      old_worklist.swap(new_worklist);
      new_worklist.clear();
    }
    if (!progress) {
      std::cerr << "[e] Infinite loop in arrival times analyzer" << std::endl;
    }
  }

  void compute_arrival_time_at_pin(signal_t const& f) {
    node_index_t n = ntk_.get_node(f);
    if (ntk_.is_pi(n)) {
      times_[f] = input_[ntk_.pi_index(n)];
    } else {
      double time = 0;
      ntk_.foreach_fanin(n, [&](auto const& fi, auto const ii) {
        time = std::max(time, times_[fi] + ntk_.get_max_pin_delay(f, ii));
      });
      times_[f] = time;
    }
  }

  void compute_arrival_time(node_index_t const& n) {
    ntk_.foreach_output(n,
                        [&](auto const& f) { compute_arrival_time_at_pin(f); });
  }
#pragma endregion

 private:
  Ntk& ntk_;
  network::incomplete_signal_map<double, Ntk> times_;
  network::tfo_manager<Ntk> tfo_;
  std::vector<double> input_;
  /* events */
  std::shared_ptr<typename mockturtle::network_events<Ntk>::add_event_type>
      add_event_;
  std::shared_ptr<typename mockturtle::network_events<Ntk>::modified_event_type>
      modified_event_;
};

}  // namespace trackers

}  // namespace analyzers

}  // namespace rinox