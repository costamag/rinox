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
  \file area_profiler.hpp
  \brief Analyzer of the area

  \author Andrea Costamagna
*/

#pragma once

#include "../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../databases/mapped_database.hpp"
#include "profilers_utils.hpp"

namespace rinox
{

namespace opto
{

namespace profilers
{

template<class Ntk, typename WinMngr>
class area_profiler
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
    cost_t mffc_area;
  };

public:
  area_profiler( Ntk& ntk, WinMngr& win_manager, profiler_params const& ps )
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

  template<class List_t>
  cost_t evaluate( List_t const& list, std::vector<signal_t> const& leaves, node_index_t const& nold = std::numeric_limits<node_index_t>::max() )
  {
    (void)nold;
    signal_t const f = insert( ntk_, leaves, list );
    node_index_t const n = ntk_.get_node( f );
    cost_t const cost_deref = recursive_deref( n );
    cost_t const cost_ref = recursive_ref( n );
    assert( std::abs( cost_ref - cost_deref ) < ps_.eps && "[e] Cost ref and deref should be the same" );
    if ( ntk_.fanout_size( n ) == 0 )
      ntk_.take_out_node( n );
    return cost_deref;
  }

  cost_t evaluate_rewiring( node_index_t const& n, std::vector<signal_t> const& new_children )
  {
    for ( auto const& f : new_children )
      ntk_.incr_fanout_size( ntk_.get_node( f ) );

    auto const& win_inputs = win_manager_.get_inputs();
    auto const cost = evaluate( n, win_inputs ) - ntk_.get_area( n );

    for ( auto const& f : new_children )
      ntk_.decr_fanout_size( ntk_.get_node( f ) );

    return cost;
  }

  cost_t evaluate( node_index_t const& n, std::vector<signal_t> const& children, node_index_t nold = std::numeric_limits<node_index_t>::max() )
  {
    std::vector<node_index_t> leaves( children.size() );
    std::transform( children.begin(), children.end(), leaves.begin(), [&]( auto const f ) { return ntk_.get_node( f ); } );
    cost_t cost_deref = measure_mffc_deref( n, leaves );
    cost_t cost_ref = measure_mffc_ref( n, leaves );
    assert( std::abs( cost_ref - cost_deref ) < ps_.eps && "[e] Cost ref and deref should be the same" );
    return cost_deref;
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn )
  {
    ntk_.foreach_gate( [&]( node_index_t const& n ) {
      if ( !ntk_.is_dead( n ) && !ntk_.is_constant( n ) && !ntk_.is_pi( n ) )
        fn( n );
    } );
  }

private:
  void compute_costs()
  {
    std::function<void( signal_t const& )> compute_costs_rec = [&]( signal_t const& f ) -> void {
      node_index_t n = ntk_.get_node( f );
      if ( ( ntk_.visited( n ) == ntk_.trav_id() ) || ntk_.is_pi( n ) )
        return;

      double const node_cost = recursive_deref( n );
      recursive_ref( n );
      assert( n < nodes_.size() );
      nodes_[n] = { n, node_cost };

      ntk_.set_visited( n, ntk_.trav_id() );
      ntk_.foreach_fanin( n, [&]( auto const& fi ) {
        compute_costs_rec( fi );
      } );
    };

    ntk_.incr_trav_id();
    ntk_.foreach_po( [&]( auto const& f ) {
      compute_costs_rec( f );
    } );
  }

  cost_t recursive_deref( node_index_t const& n )
  {
    /* terminate? */
    if ( ntk_.is_constant( n ) || ntk_.is_pi( n ) )
      return 0.0;

    /* recursively collect nodes */
    cost_t area = ntk_.get_area( n );
    ntk_.foreach_fanin( n, [&]( auto const& fi ) {
      node_index_t const ni = ntk_.get_node( fi );
      if ( ntk_.decr_fanout_size( ni ) == 0 )
      {
        area += recursive_deref( ni );
      }
    } );
    return area;
  }

  cost_t recursive_ref( node_index_t const& n )
  {
    /* terminate? */
    if ( ntk_.is_constant( n ) || ntk_.is_pi( n ) )
      return 0.0;

    /* recursively collect nodes */
    cost_t area = ntk_.get_area( n );

    ntk_.foreach_fanin( n, [&]( auto const& fi ) {
      node_index_t const ni = ntk_.get_node( fi );
      if ( ntk_.incr_fanout_size( ni ) == 0 )
      {
        area += recursive_ref( ni );
      }
    } );
    return area;
  }

  cost_t measure_mffc_deref( node_index_t const& n, std::vector<node_index_t> const& leaves )
  {
    /* reference cut leaves */
    for ( auto l : leaves )
    {
      if ( l < std::numeric_limits<uint32_t>::max() )
        ntk_.incr_fanout_size( l );
    }

    cost_t mffc_cost = recursive_deref( n );

    /* dereference leaves */
    for ( auto l : leaves )
    {
      if ( l < std::numeric_limits<uint32_t>::max() )
        ntk_.decr_fanout_size( l );
    }

    return mffc_cost;
  }

  cost_t measure_mffc_ref( node_index_t const& n, std::vector<node_index_t> const& leaves )
  {
    /* reference cut leaves */
    for ( auto l : leaves )
    {
      if ( l < std::numeric_limits<uint32_t>::max() )
        ntk_.incr_fanout_size( l );
    }

    cost_t mffc_cost = recursive_ref( n );

    /* dereference leaves */
    for ( auto l : leaves )
    {
      if ( l < std::numeric_limits<uint32_t>::max() )
        ntk_.decr_fanout_size( l );
    }

    return mffc_cost;
  }

  void sort_nodes()
  {
    nodes_.resize( ntk_.size() );
    compute_costs();
    std::stable_sort( nodes_.begin(), nodes_.end(), [&]( auto const& a, auto const& b ) {
      return a.mffc_area > b.mffc_area;
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