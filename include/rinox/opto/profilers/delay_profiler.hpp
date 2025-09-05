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
  \file delay_profiler.hpp
  \brief Analyzer of the delay

  \author Andrea Costamagna
*/

#pragma once

#include "../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../analyzers/trackers/required_times_tracker.hpp"
#include "../../evaluation/chains/bound_chain.hpp"
#include "profilers_utils.hpp"

namespace rinox
{

namespace opto
{

namespace profilers
{

template<class Ntk, typename WinMngr>
class delay_profiler
{
public:
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;
  using cost_t = double;
  static cost_t constexpr min_cost = std::numeric_limits<cost_t>::min();
  static cost_t constexpr max_cost = std::numeric_limits<cost_t>::max();
  static bool constexpr pass_window = false;
  static bool constexpr has_arrival = true;

  struct node_with_cost_t
  {
    node_index_t root;
    cost_t mffc_delay;
  };

public:
  delay_profiler( Ntk& ntk, WinMngr& win_manager, profiler_params const& ps )
      : ntk_( ntk ),
        ps_( ps ),
        nodes_( ntk_.size() ),
        arrival_( ntk_ ),
        win_manager_( win_manager )
  {
  }

  void init()
  {}

  double get_arrival( signal_t const& f ) const
  {
    return arrival_.get_time( f );
  }

  template<class List_t, class Signals>
  cost_t evaluate( List_t const& list, Signals const& leaves, node_index_t const& nold = std::numeric_limits<node_index_t>::max() )
  {
    (void)nold;
    signal_t const f = rinox::evaluation::chains::insert( ntk_, leaves, list );
    node_index_t const n = ntk_.get_node( f );
    cost_t time = 0.0;
    ntk_.foreach_output( n, [&]( auto const& f ) {
      time = std::max( time, arrival_.get_time( f ) );
    } );
    if ( ntk_.fanout_size( n ) == 0 )
      ntk_.take_out_node( n );
    return time;
  }

  template<typename Signals>
  cost_t evaluate_rewiring( node_index_t const& n, Signals const& new_children )
  {
    cost_t curr_cost = 0.0;
    ntk_.foreach_output( n, [&]( auto const& f ) {
      curr_cost = std::max( curr_cost, arrival_.get_time( f ) );
    } );
    cost_t cand_cost = 0.0;
    ntk_.foreach_output( n, [&]( auto const& f ) {
      ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
        (void)fi;
        cand_cost = std::max( cand_cost, arrival_.get_time( new_children[ii] ) + ntk_.get_max_pin_delay( f, ii ) );
      } );
    } );

    return curr_cost - cand_cost;
  }

  template<typename Signals>
  cost_t evaluate( node_index_t const& n, Signals const& children, node_index_t nold = std::numeric_limits<node_index_t>::max() )
  {
    cost_t time = 0.0;
    ntk_.foreach_output( n, [&]( auto const& f ) {
      time = std::max( time, arrival_.get_time( f ) );
    } );
    return time;
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn )
  {
    if ( ps_.max_num_roots < std::numeric_limits<uint32_t>::max() )
    {
      sort_nodes();
      uint32_t const num_roots = std::min( ps_.max_num_roots, static_cast<uint32_t>( nodes_.size() ) );

      for ( uint32_t i = 0; i < num_roots; ++i )
      {
        node_index_t const n = nodes_[i].root;
        if ( !ntk_.is_dead( n ) && !ntk_.is_constant( n ) && !ntk_.is_pi( n ) )
          fn( n );
      }
    }
    else
    {
      ntk_.foreach_gate( [&]( node_index_t const& n ) {
        if ( !ntk_.is_dead( n ) && !ntk_.is_constant( n ) && !ntk_.is_pi( n ) )
          fn( n );
      } );
    }
  }

private:
  void compute_costs()
  {
    analyzers::trackers::required_times_tracker<Ntk> required( ntk_, arrival_.worst_delay() );
    ntk_.foreach_gate( [&]( auto const& n ) {
      double node_cost = 0;
      ntk_.foreach_output( n, [&]( auto const& f ) {
        node_cost = std::max( node_cost, required.get_time( f ) - arrival_.get_time( f ) );
      } );
      assert( n < nodes_.size() );
      nodes_[n] = { n, node_cost };
    } );
  }

  void sort_nodes()
  {
    nodes_.resize( ntk_.size() );
    compute_costs();
    std::stable_sort( nodes_.begin(), nodes_.end(), [&]( auto const& a, auto const& b ) {
      return a.mffc_delay > b.mffc_delay;
    } );
  }

private:
  Ntk& ntk_;
  profiler_params const& ps_;
  std::vector<node_with_cost_t> nodes_;
  analyzers::trackers::arrival_times_tracker<Ntk> arrival_;
  WinMngr& win_manager_;
};

} /* namespace profilers */

} /* namespace opto */

} /* namespace rinox */