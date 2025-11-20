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

#ifndef RINOX_OPTO_ARRAY_AREA_MERGE_HPP_
#define RINOX_OPTO_ARRAY_AREA_MERGE_HPP_

#include <fmt/format.h>

#include <chrono>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <optional>

#include "../../../analyzers/evaluators/levels_evaluator.hpp"
#include "../../../boolean/spfd.hpp"
#include "../../../synthesis/lut/synth.hpp"
#include "../../../windowing/window_manager.hpp"
#include "../../../windowing/window_simulator.hpp"

namespace rinox {
namespace opto {
namespace array {

struct MergeAreaParamsDefault {
  static constexpr uint32_t kMaxInpsSolo = 6u;
  static constexpr uint32_t kMaxInpsDual = 5u;
  bool preserve_time = false;
};

struct MergeAreaStatsDefault {
  mockturtle::stopwatch<>::duration time{0};
  uint32_t num_matches{0};
  uint32_t num_merges{0};

  double GetTimeTotal() const { return mockturtle::to_seconds(time); }

  void RecordMerge() { num_merges++; }
  void RecordMatch() { num_matches++; }

  uint32_t GetNumMerges() const { return num_merges; }
  uint32_t GetNumMatches() const { return num_matches; }

  void Report() const {
    std::cout << fmt::format("================ TIME =====================\n");
    std::cout << fmt::format("    time total                    = {:5.2f}\n",
                             GetTimeTotal());
    std::cout << fmt::format("================ COUNT ====================\n");
    std::cout << fmt::format("    num matches   = {:5d}\n", GetNumMatches());
    std::cout << fmt::format("    num merges   = {:5d}\n", GetNumMerges());
  }
};

namespace detail {

template <typename Node>
class MatchingGraph {
  struct Edge {
    int to;
    int weight = -1;

    void Print() { std::cout << "(" << to << ", " << weight << ") "; }
  };

  struct Vertex {
    Node name;
    std::vector<Edge> edges;

    void Print() {
      std::cout << name << " : ";
      for (auto& e : edges) e.Print();
    }

    void Sort() {
      std::sort(edges.begin(), edges.end(), [&](auto const& a, auto const& b) {
        return a.weight > b.weight;
      });
    }
  };

 public:
  void AddEdge(Node i, Node j, int w) {
    if (map_.find(i) == map_.end()) {
      map_[i] = vertices_.size();
      Vertex vi{i, {}};
      vertices_.push_back(vi);
      matches_.push_back(-1);
    }

    if (map_.find(j) == map_.end()) {
      map_[j] = vertices_.size();
      Vertex vj{j, {}};
      vertices_.push_back(vj);
      matches_.push_back(-1);
    }

    Edge eij{map_[j], w};
    Edge eji{map_[i], w};

    vertices_[map_[i]].edges.push_back(eij);
    vertices_[map_[j]].edges.push_back(eji);
  }

  void Print() {
    for (int i = 0; i < vertices_.size(); ++i) {
      std::cout << i << " | ";
      vertices_[i].Print();
      std::cout << "-> " << matches_[i] << std::endl;
    }
  }

  template <typename Fn>
  void ForeachVertex(Fn&& fn) {
    for (auto& v : vertices_) v.Sort();
    std::vector<int> order(vertices_.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](auto const& a, auto const& b) {
      return vertices_[a].edges.size() < vertices_[b].edges.size();
    });
    // greedy assignment
    for (auto const& i : order) {
      if (matches_[i] >= 0) continue;
      auto const& u = vertices_[i];
      for (auto const& e : u.edges) {
        auto const& v = vertices_[e.to];
        if (matches_[e.to] < 0 && fn(u.name, v.name)) {
          matches_[e.to] = i;
          matches_[i] = e.to;
          break;
        }
      }
    }
  }

 private:
  std::unordered_map<Node, int> map_;
  std::vector<Vertex> vertices_;
  std::vector<int> matches_;
};

template <class Ntk, typename Params, typename Stats>
class MergeAreaImpl {
  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  static constexpr uint32_t kMaxNumOutputs = Ntk::kMaxNumOutputs;

 public:
  MergeAreaImpl(Ntk& ntk, Params ps, Stats& st) : ntk_(ntk), ps_(ps), st_(st) {
    for (int n = 0; n <= Params::kMaxInpsSolo; ++n) {
      proj_fns_.emplace_back();
      for (int i = 0; i < n; ++i) {
        proj_fns_[n].emplace_back(n);
        kitty::create_nth_var(proj_fns_[n][i], i);
      }
    }
  }

  void Run();
  template <typename LevelsEvaluator>
  bool ExistsPath(LevelsEvaluator&, Node, Node);
  MatchingGraph<Node> BuildMatchingGraph();
  void Merge(Node, Node);
  void MergeChild(Node, Node);

  template <typename LevelsEvaluator>
  void ChildrenCannibalism(LevelsEvaluator&);

  template <typename LevelsEvaluator>
  bool NoSidePath(Node, Signal, LevelsEvaluator&);

 private:
  Ntk& ntk_;
  Params ps_;
  Stats& st_;
  std::vector<std::vector<kitty::dynamic_truth_table>> proj_fns_;
};

// Implementation of Run split into steps
template <class Ntk, typename Params, typename Stats>
void MergeAreaImpl<Ntk, Params, Stats>::Run() {
  mockturtle::call_with_stopwatch(st_.time, [&]() {
    rinox::analyzers::trackers::arrival_times_tracker<Ntk> times(ntk_);

    auto graph = BuildMatchingGraph();
    graph.ForeachVertex([&](Node n_father, Node n_son) {
      if (ExistsPath(times, n_father, n_son)) return false;
      auto f_father = ntk_.make_signal(n_father);
      auto f_son = ntk_.make_signal(n_son);

      if (times.GetArrival(f_son) > times.GetArrival(f_father)) {
        std::swap(n_father, n_son);
        std::swap(f_father, f_son);
      }

      bool father_and_son = false;
      ntk_.foreach_fanin(n_father, [&](auto const& fi, auto ii) {
        father_and_son |= (fi == f_son);
      });
      if (father_and_son && !NoSidePath(n_father, f_son, times)) return false;

      MergeChild(n_father, n_son);
      return true;
    });
  });
}

template <class Ntk, typename Params, typename Stats>
auto MergeAreaImpl<Ntk, Params, Stats>::BuildMatchingGraph()
    -> MatchingGraph<Node> {
  std::function<void(Node)> MarkTfi = [&](Node n) {
    if (ntk_.visited(n) == ntk_.trav_id()) return;
    ntk_.set_visited(n, ntk_.trav_id());
    ntk_.foreach_fanin(
        n, [&](auto const& fi, auto) { MarkTfi(ntk_.get_node(fi)); });
  };

  MatchingGraph<Node> graph;

  analyzers::evaluators::LevelsEvaluator times(ntk_);
  std::vector<std::vector<Node>> levels;
  ntk_.foreach_gate([&](auto const& n) {
    if (!ntk_.is_dead(n) && (ntk_.num_outputs(n) == 1)) {
      auto const& f = ntk_.make_signal(n);
      int const level = times.GetArrival(f);
      while (levels.size() <= level) levels.emplace_back();
      levels[level].push_back(n);
    }
    return true;
  });

  for (int l1 = 1; l1 < levels.size(); ++l1) {
    for (auto const& n1 : levels[l1]) {
      auto const f1 = ntk_.make_signal(n1);
      ntk_.incr_trav_id();
      MarkTfi(n1);
      // start from nodes not in the TFI
      for (int l2 = l1; l2 > 0; --l2) {
        for (auto const& n2 : levels[l2]) {
          auto const& f2 = ntk_.make_signal(n2);
          if (ps_.preserve_time && times.GetArrival(f1) > times.GetRequired(f2))
            continue;
          if (ps_.preserve_time && times.GetArrival(f2) > times.GetRequired(f1))
            continue;
          if ((ntk_.visited(n2) != ntk_.trav_id())) {
            std::set<Signal> fanins;
            ntk_.foreach_fanin(n1, [&](auto const& fi, auto) {
              if (fi != f2) fanins.insert(fi);
            });
            ntk_.foreach_fanin(
                n2, [&](auto const& fi, auto) { fanins.insert(fi); });
            if (fanins.size() <= Params::kMaxInpsDual) {
              int const tA =
                  std::max(times.GetArrival(f2), times.GetArrival(f1));
              int const tR =
                  std::min(times.GetRequired(f1), times.GetRequired(f2));
              int const w = (levels.size() - tA) * (tR - tA) * fanins.size();
              graph.AddEdge(n1, n2, w);
            }
          }
        }
      }
      // nodes in the fanin
      ntk_.foreach_fanin(n1, [&](auto const& f2, auto i2) {
        auto const n2 = ntk_.get_node(f2);
        if (ntk_.num_outputs(n2) > 1) return;
        if (ntk_.is_pi(n2)) return;
        int tR = times.GetRequired(f1);
        int tA = times.GetArrival(f2);
        tA = std::max(tA, times.ArrivalWithoutFanin(f1, f2));
        tR = std::min(tR, times.RequiredWithoutFanout(f2, n1));
        if (!ps_.preserve_time || (tA <= tR)) {
          std::set<Signal> fanins;
          ntk_.foreach_fanin(n1, [&](auto const& fi, auto) {
            if (fi != f2) fanins.insert(fi);
          });
          ntk_.foreach_fanin(n2,
                             [&](auto const& fi, auto) { fanins.insert(fi); });
          if ((fanins.size() <= Params::kMaxInpsDual) &&
              NoSidePath(n1, f2, times)) {
            int const w = (levels.size() - tA) * (tR - tA) * fanins.size();
            graph.AddEdge(n1, n2, w);
          }
        }
      });
    }
  }
  return graph;
}

template <class Ntk, typename Params, typename Stats>
template <typename LevelsEvaluator>
bool MergeAreaImpl<Ntk, Params, Stats>::ExistsPath(LevelsEvaluator& times,
                                                   Node n1, Node n2) {
  if (n1 == n2) return true;
  Node n_upper = n1;
  Node n_lower = n2;
  Signal f_upper = ntk_.make_signal(n_upper);
  Signal f_lower = ntk_.make_signal(n_lower);
  double t_upper = times.GetArrival(f_upper);
  double t_lower = times.GetArrival(f_lower);

  if (t_upper < t_lower) {
    std::swap(n_upper, n_lower);
    std::swap(f_upper, f_lower);
  }

  bool father_and_son = false;
  ntk_.foreach_fanin(n_upper, [&](auto const& fi, auto ii) {
    father_and_son |= (fi == f_lower);
  });
  if (father_and_son) return false;

  ntk_.incr_trav_id();
  bool found = false;
  std::function<void(Node)> DepthFirstSearch = [&](Node u) {
    if (found) return;
    if (u == n_lower) {
      found = true;
      return;
    }
    ntk_.foreach_fanin(u, [&](auto fi, auto) {
      auto v = ntk_.get_node(fi);
      if (ntk_.visited(v) != ntk_.trav_id() &&
          (times.GetArrival(fi) >= times.GetArrival(f_lower))) {
        ntk_.set_visited(v, ntk_.trav_id());
        DepthFirstSearch(v);
      }
    });
  };
  DepthFirstSearch(n_upper);
  return found;
}

template <class Ntk, typename Params, typename Stats>
void MergeAreaImpl<Ntk, Params, Stats>::Merge(Node n1, Node n2) {
  if (ntk_.is_pi(n1) || ntk_.is_pi(n2)) return;

  auto const f2 = ntk_.make_signal(n2);
  std::unordered_map<uint64_t, std::pair<int, int>> map;
  std::set<Signal> inps;
  ntk_.foreach_fanin(n1, [&](auto const& fi, auto ii) {
    if (fi != f2) inps.insert(fi);
  });
  ntk_.foreach_fanin(n2, [&](auto const& fi, auto ii) { inps.insert(fi); });

  if (inps.size() > Params::kMaxInpsDual) return;
  int i = 0;
  std::vector<Signal> children;
  for (auto const& f : inps) {
    children.push_back(f);
    map[f.data] = std::pair{inps.size(), i++};
  }
  std::vector<kitty::dynamic_truth_table const*> sim_ptrs1;
  std::vector<kitty::dynamic_truth_table const*> sim_ptrs2;

  ntk_.foreach_fanin(n1, [&](auto const& fi, auto ii) {
    auto const p = map[fi.data];
    sim_ptrs1.push_back(&proj_fns_[p.first][p.second]);
  });
  auto funcs1 = ntk_.compute(n1, sim_ptrs1);

  ntk_.foreach_fanin(n2, [&](auto const& fi, auto ii) {
    auto const p = map[fi.data];
    sim_ptrs2.push_back(&proj_fns_[p.first][p.second]);
  });
  auto funcs2 = ntk_.compute(n2, sim_ptrs2);

  std::vector<kitty::dynamic_truth_table> tts = {funcs1[0], funcs2[0]};

  auto fnew = ntk_.create_node(children, tts);

  ntk_.substitute_node(n1, Signal{fnew.index, 0});
  ntk_.substitute_node(n2, Signal{fnew.index, 1});
}

template <class Ntk, typename Params, typename Stats>
void MergeAreaImpl<Ntk, Params, Stats>::MergeChild(Node n_father, Node n_son) {
  if (ntk_.is_pi(n_father) || ntk_.is_pi(n_son)) return;

  auto const f_son = ntk_.make_signal(n_son);
  std::unordered_map<uint64_t, std::pair<int, int>> map;
  std::set<Signal> inps;
  ntk_.foreach_fanin(n_father, [&](auto const& fi, auto ii) {
    if (fi != f_son) inps.insert(fi);
  });
  ntk_.foreach_fanin(n_son, [&](auto const& fi, auto ii) { inps.insert(fi); });

  if (inps.size() > Params::kMaxInpsDual) return;
  int i = 0;
  std::vector<Signal> children;
  for (auto const& f : inps) {
    children.push_back(f);
    map[f.data] = std::pair{inps.size(), i++};
  }
  std::vector<kitty::dynamic_truth_table const*> sim_ptrs1;
  std::vector<kitty::dynamic_truth_table const*> sim_ptrs2;

  ntk_.foreach_fanin(n_son, [&](auto const& fi, auto ii) {
    auto const p = map[fi.data];
    sim_ptrs2.push_back(&proj_fns_[p.first][p.second]);
  });
  auto funcs2 = ntk_.compute(n_son, sim_ptrs2);

  ntk_.foreach_fanin(n_father, [&](auto const& fi, auto ii) {
    auto const p = map[fi.data];
    if (fi == f_son)
      sim_ptrs1.push_back(&funcs2[0]);
    else
      sim_ptrs1.push_back(&proj_fns_[p.first][p.second]);
  });
  auto funcs1 = ntk_.compute(n_father, sim_ptrs1);

  std::vector<kitty::dynamic_truth_table> tts = {funcs1[0], funcs2[0]};

  auto fnew = ntk_.create_node(children, tts);

  ntk_.substitute_node(n_father, Signal{fnew.index, 0});
  ntk_.substitute_node(n_son, Signal{fnew.index, 1});
}

template <class Ntk, typename Params, typename Stats>
template <typename LevelsEvaluator>
bool MergeAreaImpl<Ntk, Params, Stats>::NoSidePath(Node n, Signal fi,
                                                   LevelsEvaluator& times) {
  auto min_lev = times.GetArrival(fi);
  ntk_.incr_trav_id();
  std::function<void(Signal)> Mark = [&](Signal f) {
    if (times.GetArrival(f) < min_lev) return;
    if (ntk_.is_pi(ntk_.get_node(f))) return;
    if (ntk_.visited(ntk_.get_node(f)) == ntk_.trav_id()) return;
    ntk_.set_visited(ntk_.get_node(f), ntk_.trav_id());

    ntk_.foreach_fanin(f.index, [&](auto const& fii, auto) { Mark(fii); });
  };

  ntk_.foreach_fanin(n, [&](auto const& fii, auto) {
    if (fii != fi) Mark(fii);
  });
  return (ntk_.visited(fi.index) != ntk_.trav_id());
}

template <class Ntk, typename Params, typename Stats>
template <typename LevelsEvaluator>
void MergeAreaImpl<Ntk, Params, Stats>::ChildrenCannibalism(
    LevelsEvaluator& times) {
  std::vector<std::pair<Node, Node>> merges;

  ntk_.foreach_gate([&](auto const& n) {
    if (ntk_.num_outputs(n) > 1) return true;

    int best = 0;
    std::vector<Signal> children;
    children.reserve(Params::kMaxInpsSolo);

    ntk_.foreach_fanin(
        n, [&](auto const& fi, auto /*ii*/) { children.push_back(fi); });

    for (auto const& fi : children) {
      auto const ni = ntk_.get_node(fi);
      if (ntk_.is_pi(ni)) continue;
      if (ntk_.num_outputs(ni) > 1) continue;

      if (NoSidePath(n, fi, times)) {
        std::set<Signal> fanins;

        ntk_.foreach_fanin(n, [&](auto const& fii, auto) {
          if (fii != fi) fanins.insert(fii);
        });

        ntk_.foreach_fanin(ntk_.get_node(fi),
                           [&](auto const& fii, auto) { fanins.insert(fii); });

        if (fanins.size() <= Params::kMaxInpsDual) {
          if (static_cast<int>(fanins.size()) > best) {
            best = static_cast<int>(fanins.size());
            merges.emplace_back(n, ntk_.get_node(fi));
          }
        }
      }
    }
    return true;
  });

  for (auto const& p : merges) {
    if (!ntk_.is_dead(p.first) && !ntk_.is_dead(p.second))
      MergeChild(p.first, p.second);
  }
}

}  // namespace detail

// Public API
template <class Ntk, typename Params = MergeAreaParamsDefault,
          typename Stats = MergeAreaStatsDefault>
void MergeArea(Ntk& ntk, Params ps = {}, Stats* pst = nullptr) {
  Stats st;
  detail::MergeAreaImpl<Ntk, Params, Stats> p(ntk, ps, st);
  p.Run();
  if (pst) *pst = st;
}

}  // namespace array
}  // namespace opto
}  // namespace rinox

#endif  // RINOX_OPTO_ARRAY_AREA_MERGE_HPP_