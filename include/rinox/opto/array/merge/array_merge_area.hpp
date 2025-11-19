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

#ifndef RINOX_ALGORITHMS_ARRAY_MERGE_AREA_HPP_
#define RINOX_ALGORITHMS_ARRAY_MERGE_AREA_HPP_

#include <fmt/format.h>

#include <chrono>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <optional>

#include "../../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../../analyzers/evaluators/levels_evaluator.hpp"
#include "../../../boolean/simd.hpp"
#include "../../../boolean/spfd.hpp"
#include "../../../dependency/window/oracles.hpp"
#include "../../../network/visualize.hpp"
#include "../../../synthesis/lut/synth.hpp"
#include "../../../windowing/simula_manager.hpp"
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
namespace array {

struct MergeAreaParamsDefault {
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr uint32_t kMaxInpsSolo = 6u;
  bool preserve_time = false;
};

struct MergeAreaStatsDefault {
  mockturtle::stopwatch<>::duration time_matching{0};
  mockturtle::stopwatch<>::duration time_covering{0};
  uint32_t num_merges{0};

  double GetTimeMatching() const {return mockturtle::to_seconds(time_matching);}
  double GetTimeCovering() const {return mockturtle::to_seconds(time_covering);}

  double GetTimeTotal() const { return  GetTimeMatching() +
                                        GetTimeCovering(); }

  void RecordMerge() { num_merges++; }

  uint32_t GetNumMerges() const { return num_merges; }

  void Report() const {
    std::cout << fmt::format("================ TIME  ====================\n");
    std::cout << fmt::format("    time matching                 = {:5.2f}\n",
                             GetTimeMatching());
    std::cout << fmt::format("    time covering                 = {:5.2f}\n",
                             GetTimeCovering());
    std::cout << fmt::format("    time total                    = {:5.2f}\n",
                             GetTimeTotal());
    std::cout << fmt::format("================ QUALITY ==================\n");
    std::cout << fmt::format("    count merges                  = {:5d}\n",
                             GetNumMerges());
  }
};

namespace detail {

template<typename Node>
class MatchingGraph{

    struct Edge{
        int to;
        int weight=-1;

        void Print()
        {
            std::cout << "(" << to << ", " << weight << ") ";
        }
    };

    struct Vertex{
      Node name;
      std::vector<Edge> edges;

      void Print()
      {
        std::cout << name << " : ";
        for ( auto & e : edges )
            e.Print();
      }

      void Sort()
      {
        std::sort( edges.begin(), edges.end(), [&](auto const& a, auto const& b){
            return a.weight > b.weight;
        } );
      }
    };

public:

    void Clear()
    {
      map_.clear();
      vertices_.clear();
      matches_.clear();
    }

    void AddEdge(Node i, Node j, int w)
    {
        if ( map_.find(i) == map_.end() )
        {
            map_[i]=vertices_.size();
            Vertex vi{i, {}};
            vertices_.push_back(vi);
            matches_.push_back(-1);
        }

        if ( map_.find(j) == map_.end() )
        {
            map_[j]=vertices_.size();
            Vertex vj{j, {}};
            vertices_.push_back(vj);
            matches_.push_back(-1);
        }

        Edge eij{ map_[j], w };
        Edge eji{ map_[i], w };

        vertices_[map_[i]].edges.push_back(eij);
        vertices_[map_[j]].edges.push_back(eji);
    }

    void Print()
    {
        for ( int i=0; i<vertices_.size(); ++i )
        {
            std::cout << i << " | ";
            vertices_[i].Print();
            std::cout << "-> " << matches_[i] << std::endl;
        }
    }


    template<typename Fn>
    void ForeachMatch( Fn && fn )
    {
      for (auto& v: vertices_)
          v.Sort();
      std::vector<int> order(vertices_.size());
      std::iota(order.begin(), order.end(),0);
      std::sort(order.begin(), order.end(), [&](auto const& a, auto const& b){
          return vertices_[a].edges.size() < vertices_[b].edges.size();
      });
      // greedy assignment
      for (auto const& i : order)
      {
        if(matches_[i] >= 0) continue;
        auto const& u = vertices_[i];
        for (auto const& e : u.edges)
        {
          auto const& v = vertices_[e.to];
          if(matches_[e.to] < 0 && fn(u.name, v.name))
          {
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
  MergeAreaImpl(Ntk& ntk, Params ps, Stats& st)
      : // parameters
        ps_(ps),
        // statistics
        st_(st),
        // attributes
        ntk_(ntk),
        levels_(ntk)
  {
    for (int n = 0; n <= Params::kMaxInpsSolo; ++n) {
      proj_fns_.emplace_back();
      for (int i = 0; i < n; ++i) {
        proj_fns_[n].emplace_back(n);
        kitty::create_nth_var(proj_fns_[n][i], i);
      }
    }
  }

  void Report();
  void Run();

 private:
  bool SkipNode(Node) const;
  void Matching(int);
  bool Covering();
  void MergeChild(Node, Node);
  bool NoSidePath(Node, Signal);
  bool ExistsPath(Node, Node);

 private:
  // Parameters
  Params ps_;
  // Statistics
  Stats& st_;
 // Attributes
  Ntk& ntk_;
  rinox::analyzers::evaluators::LevelsEvaluator<Ntk> levels_;
  MatchingGraph<Node> graph_;
  std::vector<std::vector<kitty::dynamic_truth_table>> proj_fns_;

};

template <typename Ntk, typename Params, typename Stats>
void MergeAreaImpl<Ntk, Params, Stats>::Report() {
}

// Implementation of Run split into steps
template <typename Ntk, typename Params, typename Stats>
void MergeAreaImpl<Ntk, Params, Stats>::Run() {

  int const last_level = levels_.GetLast();

  int level = 1;
  while(level<=last_level)
  {
    mockturtle::call_with_stopwatch(st_.time_matching, [&]() { Matching(level); });
    bool const update = mockturtle::call_with_stopwatch(st_.time_covering, [&]() { return Covering(); });
    if(update) levels_.ComputeLevels();
    ++level;
  }
}

template <typename Ntk, typename Params, typename Stats>
void MergeAreaImpl<Ntk, Params, Stats>::Matching(int level)
{
  std::function<void(Node)> MarkTfi = [&](Node n){
    if(ntk_.visited(n)==ntk_.trav_id()) return;
    ntk_.set_visited(n, ntk_.trav_id());
    ntk_.foreach_fanin(n, [&](auto const& fi, auto){
      MarkTfi(ntk_.get_node(fi));
    });
  };

  graph_.Clear();
  levels_.ForeachNode(level, [&](Node n1){
    if(ntk_.num_outputs(n1)>1) return true;
    auto const f1 = ntk_.make_signal(n1);
    ntk_.incr_trav_id();
    MarkTfi(n1);
    // start from nodes not in the TFI
    for(int l2=levels_.GetLast(); l2>0; --l2)
    {
      levels_.ForeachNode(level, [&](Node n2){
        if(ntk_.num_outputs(n2)>1) return true;
        auto const& f2 = ntk_.make_signal(n2);
        if(ps_.preserve_time && levels_.GetArrival(f1)>levels_.GetRequired(f2)) return true;
        if(ps_.preserve_time && levels_.GetArrival(f2)>levels_.GetRequired(f1)) return true;
        if((ntk_.visited(n2)!=ntk_.trav_id()))
        {
            std::set<Signal> fanins;
            ntk_.foreach_fanin(n1, [&](auto const& fi, auto){ if(fi!=f2) fanins.insert(fi);});
            ntk_.foreach_fanin(n2, [&](auto const& fi, auto){ fanins.insert(fi);});
            if (fanins.size() <= Params::kMaxInpsDual)
            {
              int const tA = std::max(levels_.GetArrival(f2), levels_.GetArrival(f1));
              int const tR = std::min(levels_.GetRequired(f1), levels_.GetRequired(f2));
              int const w = (level-tA)*(tR-tA)*fanins.size();
              graph_.AddEdge( n1, n2, w );
            }
        }
        return true;
      });
    }
    // nodes in the fanin
    ntk_.foreach_fanin(n1, [&](auto const& f2, auto i2){
        auto const n2 = ntk_.get_node(f2);
        if(ntk_.num_outputs(n2)>1) return;
        if(ntk_.is_pi(n2)) return;
        int tR = levels_.GetRequired(f1);
        int tA = levels_.GetArrival(f2);
        tA = std::max(tA, levels_.ArrivalWithoutFanin(f1, f2));
        tR = std::min(tR, levels_.RequiredWithoutFanout(f2, n1));
        if(!ps_.preserve_time || (tA<=tR))
        {
          std::set<Signal> fanins;
          ntk_.foreach_fanin(n1, [&](auto const& fi, auto){ if(fi!=f2) fanins.insert(fi);});
          ntk_.foreach_fanin(n2, [&](auto const& fi, auto){ fanins.insert(fi);});
          if ((fanins.size() <= Params::kMaxInpsDual) && NoSidePath(n1, f2))
          {
            int const w = (level-tA)*(tR-tA)*fanins.size();
            graph_.AddEdge( n1, n2, w );
          }
        }
    });
    return true;
  });

}

template <typename Ntk, typename Params, typename Stats>
bool MergeAreaImpl<Ntk, Params, Stats>::Covering()
{
    bool update=false;
    graph_.ForeachMatch( [&](Node n_father, Node n_son ){
      if (ExistsPath(n_father, n_son)) return false;
      auto f_father = ntk_.make_signal(n_father);
      auto f_son = ntk_.make_signal(n_son);
        
      if(levels_.GetArrival(f_son)>levels_.GetArrival(f_father))
      {
        std::swap(n_father, n_son);
        std::swap(f_father, f_son);
      }
    
      bool father_and_son=false;
      ntk_.foreach_fanin(n_father, [&](auto const& fi, auto ii){
        father_and_son |= (fi==f_son);
      });
      if(father_and_son && !NoSidePath(n_father, f_son)) return false;
    
      MergeChild(n_father, n_son);
      update=true;
      return true;
    } );
    return update;
}

template <typename Ntk, typename Params, typename Stats>
void MergeAreaImpl<Ntk, Params, Stats>::MergeChild(Node n_father, Node n_son) {
  if (ntk_.is_pi(n_father) || ntk_.is_pi(n_son)) return;

  auto const f_son = ntk_.make_signal(n_son);
  std::unordered_map<uint64_t, std::pair<int, int>> map;
  std::set<Signal> inps;
  ntk_.foreach_fanin(n_father, [&](auto const& fi, auto ii) { if(fi!=f_son) inps.insert(fi);});
  ntk_.foreach_fanin(n_son, [&](auto const& fi, auto ii) { inps.insert(fi);});

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
    if(fi==f_son)
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

template <typename Ntk, typename Params, typename Stats>
bool MergeAreaImpl<Ntk, Params, Stats>::SkipNode(Node n) const {
  if (ntk_.fanout_size(n) <= 0) return true;
  if (ntk_.is_pi(n) || ntk_.is_constant(n)) return true;
  if (ntk_.is_dead(n)) return true;
  if (ntk_.num_outputs(n)>1) return true;
  return false;
}

template <typename Ntk, typename Params, typename Stats>
bool MergeAreaImpl<Ntk, Params, Stats>::NoSidePath(Node n, Signal fi)
{
  auto min_lev = levels_.GetArrival(fi);
  ntk_.incr_trav_id();
  std::function<void(Signal)> Mark=[&](Signal f){
    if(levels_.GetArrival(f)<min_lev) return;
    if(ntk_.is_pi(ntk_.get_node(f))) return;
    if(ntk_.visited(ntk_.get_node(f))==ntk_.trav_id()) return;

    ntk_.set_visited(ntk_.get_node(f), ntk_.trav_id());
    ntk_.foreach_fanin(f.index, [&](auto const & fii, auto){
      Mark(fii);
    });
  };

  ntk_.foreach_fanin(n, [&](auto const & fii, auto){
    if(fii!=fi) Mark(fii);
  });
  return(ntk_.visited(fi.index)!=ntk_.trav_id());
}

template <typename Ntk, typename Params, typename Stats>
bool MergeAreaImpl<Ntk, Params, Stats>::ExistsPath(Node n1, Node n2) {
  if (n1 == n2) return true;
  Node n_upper = n1;
  Node n_lower = n2;
  Signal f_upper = ntk_.make_signal(n_upper);
  Signal f_lower = ntk_.make_signal(n_lower);
  double t_upper = levels_.GetArrival(f_upper);
  double t_lower = levels_.GetArrival(f_lower);

  if (t_upper < t_lower) {
    std::swap(n_upper, n_lower);
    std::swap(f_upper, f_lower);
  }

  bool father_and_son=false;
  ntk_.foreach_fanin(n_upper, [&](auto const& fi, auto ii){
    father_and_son |= (fi==f_lower);
  });
  if(father_and_son) return false;

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
          (levels_.GetArrival(fi) >= levels_.GetArrival(f_lower))) {
        ntk_.set_visited(v, ntk_.trav_id());
        DepthFirstSearch(v);
      }
    });
  };
  DepthFirstSearch(n_upper);
  return found;
}

}  // namespace detail

// Public API
template <class Ntk, typename Params = MergeAreaParamsDefault, typename Stats = MergeAreaStatsDefault>
void MergeArea( Ntk& ntk,
                      Params ps = {},
                      Stats* pst = nullptr) {
  Stats st;
  detail::MergeAreaImpl<Ntk, Params, Stats> p(ntk, ps, st);
  p.Run();
  st.Report();
  if (pst) *pst = st;
}

}  // namespace array
}  // namespace opto
}  // namespace rinox

#endif  // RINOX_ALGORITHMS_ARRAY_AREA_ORCHESTRA_HPP_
