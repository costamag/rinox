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
 * \file _Levels_eval.hpp
 * \brief Compute  time information of a network and updates upon change
 *
 * \author Andrea Costamagna
 */

#pragma once

#include "../../network/signal_map.hpp"
#include "../../network/tfo_manager.hpp"
#include <limits>
#include <mockturtle/networks/events.hpp>

namespace rinox
{

namespace analyzers
{

namespace evaluators
{

/*! \brief Engine to evaluate the  Levels of a network.
 *
 * This engine computes the  Levels of a network and keeps them up-to-date.
 * During construction it is possible to specify a vector of  Levels. If
 * not provided, the engine assumes 0  time at all the PIs. At construction,
 * the  Levels are propagated in the network.
 *
 * The engine is equipped with a transitive fanout (TFO) manager to handle the
 * timing updates. Two events trigger the update:
 * - Node addition: The  time of the node is computed from the fanins
 * - Node modification: The  time of the TFO of the affected nodes is updated.
 *
 * \tparam Ntk the network type to be analyzed.
 *
 * \verbatim embed:rst

   Example

   .. code-block:: c++

      bound_network ntk( gates );
      auto const a = ntk.create_pi();
      auto const f1 = ntk.create_node( { a }, 0 );
      _Levels_eval ( ntk );
      auto const f2 = ntk.create_node( { f1 }, 0 );
   \endverbatim
 */
template<class Ntk>
class LevelsEvaluator
{
public:
  using Index = typename Ntk::node_index_t;
  using Node = typename Ntk::node_index_t;
  using Signal = typename Ntk::signal_t;

public:
  LevelsEvaluator( Ntk& ntk )
      : ntk_( ntk )
  {
    ComputeLevels();
  }

  int GetLast()
  {
    return levels_.size()-1;
  }

  template<typename Fn>
  void ForeachNode(int level, Fn && fn)
  {
    for(auto const& n : levels_[level])
    {
      if(!fn(n)) return;
    }
  }

  void ComputeLevels()
  {
    levels_.clear();
    levels_.emplace_back();
    arrivals_.clear();
    required_.clear();
    ntk_.foreach_pi([&](auto const& n){
      arrivals_[n]=0;
      levels_[0].push_back(n);
      required_[n]=0;
    });
    ComputeArrivalLevels();
    ComputeRequiredLevels();
  }
#pragma region Interface methods
public:

  [[nodiscard]] int GetArrival( Signal const f ) { return arrivals_[ntk_.get_node(f)]; }
  [[nodiscard]] int GetRequired( Signal const f ) { return required_[ntk_.get_node(f)]; }

#pragma endregion

private:
  #pragma region Arrivals

    void ComputeArrivalLevels()
    {
      ntk_.incr_trav_id();
      ntk_.foreach_pi([&](auto const& n){ ntk_.set_visited(n, ntk_.trav_id()); });
      ntk_.foreach_po([&]( auto const& f) {
        ComputeArrivalLevelsTfi(f);
      } );
    }

void ComputeArrivalLevelsTfi( Signal const& f )
{
  auto const n = ntk_.get_node( f );

  if ( ntk_.visited( n ) == ntk_.trav_id() )
    return;
  if ( ntk_.is_dead( n ) )
    return;

  // Mark as visited immediately to break cycles / repeated visits
  ntk_.set_visited( n, ntk_.trav_id() );

  ntk_.foreach_fanin( n, [&]( auto const& fi, auto ) {
    ComputeArrivalLevelsTfi( fi );
  } );

  ComputeArrivalTime( n );
}

void ComputeArrivalTime( Index const& n )
{
  arrivals_[n] = 0;
  ntk_.foreach_fanin( n, [&](auto const& fi, auto) {
    arrivals_[n] = std::max( arrivals_[n], arrivals_[ntk_.get_node(fi)] + 1 );
  } );
  // no need to set_visited here anymore
  while ( levels_.size() <= static_cast<size_t>(arrivals_[n]) )
    levels_.emplace_back();
  levels_[arrivals_[n]].push_back( n );
}


#pragma endregion

  #pragma region Required

  void ComputeRequiredLevels()
  {
    int const max_level = levels_.size()-1;
    required_.clear();

    for(int l=levels_.size()-1; l>=0; --l)
    {
      for(auto const& n : levels_[l])
      {
        required_[n]=max_level;
        ntk_.foreach_fanout(n,[&](auto const& no){
          required_[n]=std::min(required_[n], required_[no]-1);
        });
      }
    }
  }

public:
  int ArrivalWithoutFanin( Signal const& f, Signal const& child )
  {
    Index n = ntk_.get_node( f );
    int time = 0;
    ntk_.foreach_fanin( n, [&]( auto const& fi, auto const ii ) {
      if(fi!=child)
      {
        time = std::max( time, arrivals_[ntk_.get_node(fi)] + 1 );
      }
    } );
    return time;
  }

  int RequiredWithoutFanout( Signal const& f, Index const& father )
  {
    auto time = std::numeric_limits<int>::max();
    auto const n = ntk_.get_node(f);

    ntk_.foreach_fanout(n, [&](auto const& no) {
      if(no != father)
      {
        ntk_.foreach_fanin(no, [&](auto const& fi, auto ii) {
          if (fi == f)
          {
            time = std::min(time, required_[no] - 1);
          }
        });
      }
    });
    return time;
  }

#pragma endregion

private:
  Ntk& ntk_;
  std::unordered_map<Node, int> arrivals_;
  std::unordered_map<Node, int> required_;
  std::vector<std::vector<Node>> levels_;
};


} // namespace evals
} // namespace analyzers
} // namespace rinox