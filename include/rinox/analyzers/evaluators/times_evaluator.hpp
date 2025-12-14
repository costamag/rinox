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
 * \file times_evaluator.hpp
 * \brief Compute  time information of a network and updates upon change
 *
 * \author Andrea Costamagna
 */

#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <mockturtle/networks/events.hpp>
#include <unordered_map>
#include <vector>

#include "../../network/signal_map.hpp"

namespace rinox {

namespace analyzers {

namespace evaluators {

/*! \brief Engine to evaluate the levels of a network and keep them updated
 * within slices
 *
 */
template <class Ntk>
class TimesEvaluator {
 public:
  using Node = typename Ntk::node_index_t;
  using Signal = typename Ntk::signal_t;

  struct Info {
    Info(Node node) : node(node) {}
    /* actual information */
    Node node;
    /* structural */
    int32_t tier = 0;
    Info* prev = nullptr;
    Info* next = nullptr;
    friend struct InfoList;
  };

  struct InfoList {
    int32_t tier = -1;
    Info* tail = nullptr;
    Info* head = nullptr;

    InfoList(int32_t tier) : tier(tier) {}

    void Add(Info* p) {
      p->prev = head;
      p->next = nullptr;
      p->tier = tier;
      if (head) head->next = p;
      head = p;
      if (!tail) tail = p;
    }

    void Remove(Info* p) {
      if (p->tier != tier) return;
      if (p == tail && p == head) {
        tail = nullptr;
        head = nullptr;
        p->prev = nullptr;
        p->next = nullptr;
        p->tier = -1;
      } else if (p == tail) {
        tail = p->next;
        tail->prev = nullptr;
        p->next = nullptr;
        p->prev = nullptr;
        p->tier = -1;
      } else if (p == head) {
        head = p->prev;
        head->next = nullptr;
        p->next = nullptr;
        p->prev = nullptr;
        p->tier = -1;
      } else {
        p->prev->next = p->next;
        p->next->prev = p->prev;
        p->tier = -1;
      }
    }

    template <typename Fn>
    void ForeachEntry(Fn&& fn) {
      Info* p = tail;
      Info* next;
      while (p) {
        next = p->next;
        fn(p);
        p = next;
      }
    }
  };

 public:
  TimesEvaluator(Ntk& ntk) : ntk_(ntk), arrivals_(ntk), required_(ntk) {
    // input tier
    InfoList list(0);
    ntk_.foreach_pi([&](auto const& n) {
      auto info_ptr = std::make_unique<Info>(n);
      Info* p = info_ptr.get();
      list.Add(p);
      info_.emplace(n, std::move(info_ptr));
    });
    tiers_.push_back(std::move(list));
    Init();

    Run();
  }

  ~TimesEvaluator() {
    if (add_event_) {
      ntk_.events().release_add_event(add_event_);
    }
  }

  void Init() {
    add_event_ = ntk_.events().register_add_event([&](const Node& n) {
      ComputeTiers(n);
      ComputeArrival(n);
    });
  }

#pragma region Interface methods
 public:
  [[nodiscard]] int GetTier(Node const& n) {
    if (info_.find(n) == info_.end()) {
      ComputeTiers(n);
      ComputeArrival(n);
    }
    return info_[n]->tier;
  }
  [[nodiscard]] int GetTier(Signal const& f) {
    return GetTier(ntk_.get_node(f));
  }
  [[nodiscard]] double GetArrival(Signal const& f) { return arrivals_[f]; }
  [[nodiscard]] double GetRequired(Signal const& f) { return required_[f]; }
  [[nodiscard]] double GetRequired(Node const& n) {
    double required = std::numeric_limits<double>::max();
    ntk_.foreach_output(
        n, [&](auto const& f) { required = std::min(required, required_[f]); });
    return required;
  }
#pragma endregion

  void Run() {
    ComputeTiers();
    ComputeTimes();
  }

#pragma region Tiers

  void ComputeTiers() {
    ntk_.incr_trav_id();
    ntk_.foreach_pi(
        [&](auto const& n) { ntk_.set_visited(n, ntk_.trav_id()); });
    ntk_.foreach_po([&](auto const& f) { ComputeTiersTfi(f); });
  }

  void ComputeTiersTfi(Signal const& f) {
    auto const n = ntk_.get_node(f);
    if (ntk_.visited(n) == ntk_.trav_id()) return;
    if (ntk_.is_dead(n)) return;

    ntk_.set_visited(n, ntk_.trav_id());
    ntk_.foreach_fanin(n, [&](auto const& fi, auto) { ComputeTiersTfi(fi); });
    ComputeTiers(n);
  }

  void ComputeTiers(Node const& n) {
    if (ntk_.is_pi(n)) return;

    Info* p = nullptr;
    auto it = info_.find(n);
    if (it == info_.end()) {
      auto info_ptr = std::make_unique<Info>(n);
      p = info_ptr.get();
      info_.emplace(n, std::move(info_ptr));
    } else {
      p = it->second.get();
      tiers_[p->tier].Remove(p);
    }

    int32_t tier = -1;
    ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
      Node const ni = ntk_.get_node(fi);
      Info* pi = info_[ni].get();
      tier = std::max(tier, pi->tier + 1);
    });

    while (static_cast<int32_t>(tiers_.size()) <= tier)
      tiers_.emplace_back(static_cast<int32_t>(tiers_.size()));
    tiers_[tier].Add(p);
  }

#pragma endregion

  void UpdateTiers(Node n, int32_t bottom, int32_t top) {
    ntk_.incr_trav_id();
    std::vector<Node> tfo;
    CollectTfo(tfo, n, top);

    ntk_.incr_trav_id();
    std::vector<Node> tfi;
    CollectTfi(tfi, n, bottom);

    std::sort(tfo.begin(), tfo.end(),
              [&](Node a, Node b) { return info_[a]->tier > info_[b]->tier; });
    std::sort(tfi.begin(), tfi.end(),
              [&](Node a, Node b) { return info_[a]->tier > info_[b]->tier; });

    for (int i = static_cast<int>(tfo.size()) - 1; i >= 0; --i) {
      ComputeTiers(tfo[i]);
      ComputeArrival(tfo[i]);
    };
    for (int i = 0; i < tfo.size(); ++i) ComputeRequired(tfo[i]);
    for (int i = 0; i < tfi.size(); ++i) ComputeRequired(tfi[i]);
  }

  void CollectTfi(std::vector<Node>& tfi, Node n, int bottom) {
    if (ntk_.visited(n) == ntk_.trav_id()) return;
    if (info_[n]->tier < bottom) return;
    ntk_.set_visited(n, ntk_.trav_id());
    tfi.push_back(n);
    ntk_.foreach_fanin(n, [&](auto const& fi, auto i) {
      Node ni = ntk_.get_node(fi);
      CollectTfi(tfi, ni, bottom);
    });
  }

  void CollectTfo(std::vector<Node>& tfo, Node n, int top) {
    if (ntk_.visited(n) == ntk_.trav_id()) return;
    if (info_[n]->tier > top) return;
    ntk_.set_visited(n, ntk_.trav_id());
    tfo.push_back(n);
    ntk_.foreach_fanout(n, [&](auto const& no) { CollectTfo(tfo, no, top); });
  }

  void ComputeTimes() {
    ComputeArrivals();
    ComputeRequired();
  }

  template <typename Fn>
  void ForeachNode(int32_t tier, Fn&& fn) {
    auto& list = tiers_[tier];
    Info* p = list.tail;
    while (p) {
      Info* next = p->next;
      fn(p->node);
      p = next;
    }
  }

  template <typename Fn>
  void ForeachTierReverse(Fn&& fn) {
    int32_t const bottom = 0;
    int32_t const top = static_cast<int32_t>(tiers_.size()) - 1;
    ForeachTierReverse(bottom, top, fn);
  }

  template <typename Fn>
  void ForeachTierReverse(int32_t bottom, int32_t top, Fn&& fn) {
    for (int32_t tier = top; tier >= bottom; --tier) fn(tiers_[tier]);
  }

  template <typename Fn>
  void ForeachTierForward(Fn&& fn) {
    int32_t const bottom = 0;
    int32_t const top = static_cast<int32_t>(tiers_.size()) - 1;
    ForeachTierForward(bottom, top, fn);
  }

  template <typename Fn>
  void ForeachTierForward(int32_t bottom, int32_t top, Fn&& fn) {
    for (int32_t tier = bottom; tier <= top; ++tier) fn(tiers_[tier]);
  }

  void ComputeArrivals() {
    arrivals_.reset();
    ForeachTierForward([&](InfoList& node_list) {
      node_list.ForeachEntry([&](Info* p) { ComputeArrival(p->node); });
    });
  }

  void ComputeArrival(Node n) {
    ntk_.foreach_output(n, [&](auto const& f) {
      bool any_fanin = false;
      double val = std::numeric_limits<double>::lowest();
      ntk_.foreach_fanin(f, [&](auto const& fi, auto i) {
        any_fanin = true;
        val = std::max(val, arrivals_[fi] + ntk_.get_max_pin_delay(f, i));
      });
      arrivals_[f] = any_fanin ? val : 0.0;  // PIs / no-fanin outputs get 0
    });
  }

  void ComputeRequired(Node n) {
    ntk_.foreach_output(n, [&](auto const& f) {
      bool any_fanout = false;
      double val = std::numeric_limits<double>::max();
      ntk_.foreach_fanout(f, [&](auto const& no) {
        ntk_.foreach_output(no, [&](auto const& fo) {
          ntk_.foreach_fanin(fo, [&](auto const& fi, auto i) {
            if (fi == f) {
              val =
                  std::min(val, required_[fo] - ntk_.get_max_pin_delay(fo, i));
              any_fanout = true;
            }
          });
        });
      });
      required_[f] = any_fanout ? val : *worst_arrival_;
    });
  }

  double GetWorstArrival() {
    double res = 0;
    auto& list = tiers_.back();
    list.ForeachEntry([&](Info* p) {
      ntk_.foreach_output(
          p->node, [&](auto const& f) { res = std::max(res, arrivals_[f]); });
    });
    return res;
  }

  int32_t GetWorstTier() {
    int32_t res = 0;
    ntk_.foreach_po([&](auto const& f) {
      auto n = ntk_.get_node(f);
      res = std::max(res, info_[n]->tier);
    });
    return res;
  }

#pragma region Required
  void ComputeRequired() {
    required_.reset();
    double const worst_arrival = GetWorstArrival();
    if (!worst_arrival_) worst_arrival_ = worst_arrival;
    auto& list = tiers_.back();

    ntk_.incr_trav_id();
    list.ForeachEntry([&](Info* p) {
      ntk_.foreach_output(p->node,
                          [&](auto const& f) { required_[f] = worst_arrival; });
      ntk_.set_visited(p->node, ntk_.trav_id());
    });

    ForeachTierReverse([&](InfoList& node_list) {
      node_list.ForeachEntry([&](Info* p) {
        ntk_.foreach_output(p->node, [&](auto const& f) {
          ntk_.foreach_fanin(f, [&](auto const& fi, auto i) {
            Node ni = ntk_.get_node(fi);
            if (ntk_.visited(ni) != ntk_.trav_id()) {
              ntk_.set_visited(ni, ntk_.trav_id());
              required_[fi] = std::numeric_limits<double>::max();
            }
            required_[fi] = std::min(
                required_[fi], required_[f] - ntk_.get_max_pin_delay(f, i));
          });
        });
      });
    });
  }

#pragma endregion

 private:
  Ntk& ntk_;
  std::optional<double> worst_arrival_;
  rinox::network::incomplete_signal_map<double, Ntk> arrivals_;
  rinox::network::incomplete_signal_map<double, Ntk> required_;
  std::unordered_map<Node, std::unique_ptr<Info>> info_;
  std::vector<InfoList> tiers_;
  std::shared_ptr<typename mockturtle::network_events<Ntk>::add_event_type>
      add_event_;
};

}  // namespace evaluators
}  // namespace analyzers
}  // namespace rinox