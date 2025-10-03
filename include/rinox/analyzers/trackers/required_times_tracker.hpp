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
 * \file required_times_tracker.hpp
 * \brief Compute required time information of a network and updates it upon change
 *
 * \author Andrea Costamagna
 */

#pragma once

#include "../../network/signal_map.hpp"
#include "topo_sort_tracker.hpp"
#include <limits>

namespace rinox
{

namespace analyzers
{

namespace trackers
{

/*! \brief Engine to evaluate the required times of a network.
 *
 * This engine computes the required times of a network and keeps them up-to-date.
 * During construction it is possible to specify a vector of required times. At
 * construction, the required times are propagated in the network.
 *
 * The engine is equipped with an engine which maintains the topological order
 * of the network up-to-date. Required times computation can be extremely
 * expensive when called multiple times during graph optimization. The on_modified
 * event ie equipped with techniques to minimizing updates. The engine maintaining
 * the topological order serves the same purpose.
 *
 * \tparam Ntk the network type to be analyzed.
 *
 * \verbatim embed:rst

   Example

   .. code-block:: c++

      bound_network ntk( gates );
      auto const a = ntk_.create_pi();
      auto const f1 = ntk_.create_node( { a }, 0 );
      required_times_tracker required( ntk );
      auto const f2 = ntk_.create_node( { f1 }, 0 );
   \endverbatim
 */
template<class Ntk>
class required_times_tracker
{
public:
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;
  static constexpr double infinite_time = std::numeric_limits<double>::max();

public:
  required_times_tracker( Ntk& ntk, double const& required )
      : ntk_( ntk ),
        times_( ntk ),
        topo_sort_( ntk ),
        output_( std::vector<double>( ntk_.num_pos(), required ) )
  {
    init();
  }

  required_times_tracker( Ntk& ntk, std::vector<double> const& output_required )
      : ntk_( ntk ),
        times_( ntk ),
        topo_sort_( ntk ),
        output_( output_required )
  {
    init();
  }

  void init()
  {
    /* check if the network implements the needed functionalities */
    static_assert( mockturtle::is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( mockturtle::has_incr_trav_id_v<Ntk>, "Ntk does not implement the incr_trav_id method" );
    static_assert( mockturtle::has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
    static_assert( rinox::traits::has_foreach_output_v<Ntk>, "Ntk does not implement the foreach_output method" );
    static_assert( mockturtle::has_foreach_fanout_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
    static_assert( mockturtle::has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
    static_assert( mockturtle::has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_pi method" );
    static_assert( mockturtle::has_set_value_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( mockturtle::has_value_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( mockturtle::has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( mockturtle::has_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( mockturtle::has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
    static_assert( mockturtle::has_size_v<Ntk>, "Ntk does not implement the size method" );

    compute_required_times();

    add_event_ = ntk_.events().register_add_event( [&]( const node_index_t& n ) {
      times_.resize();
      ntk_.foreach_output( n, [&]( auto f ) {
        times_[f] = infinite_time;
      } );
    } );

    modified_event_ = ntk_.events().register_modified_event( [&]( const auto& n, auto old_children ) {
      /* mark the TFI by setting the required time to infinite ( unless PO ) */
      ntk_.incr_trav_id();
      if ( has_required_time_update( n ) )
      {
        reset_required_tfi( n );
      }

      auto new_children = ntk_.get_children( n );
      for ( auto const& f : new_children )
      {
        if ( has_required_time_update( f ) )
          reset_required_tfi( ntk_.get_node( f ) );
      }

      for ( auto const& f : old_children )
      {
        if ( has_required_time_update( f ) )
          reset_required_tfi( ntk_.get_node( f ) );
      }
      /* now all the nodes to be modified have the visited flag set to trav_id */

      /* find the highest level */
      auto level = topo_sort_.get_level( n );
      for ( auto const& f : old_children )
      {
        level = std::max( level, topo_sort_.get_level( ntk_.get_node( f ) ) );
      }

      /* update the affected required times from one level up */
      update_required( level );
    } );
  }

  ~required_times_tracker()
  {
    if ( add_event_ )
    {
      ntk_.events().release_add_event( add_event_ );
    }

    if ( modified_event_ )
    {
      ntk_.events().release_modified_event( modified_event_ );
    }
  }

  void SetMaxDelay( double delay )
  {
    uint32_t i = 0;
    ntk_.foreach_po( [&]( auto f ) {
      output_[i++] = delay;
    } );
    compute_required_times();
  }

#pragma region Interface methods
public:
  [[nodiscard]] double get_time( signal_t const f ) const
  {
    return times_[f];
  }
#pragma endregion

#pragma region Implementation details
private:
  [[nodiscard]] bool is_marked_todo( node_index_t const& n ) const
  {
    return ntk_.visited( n ) == ntk_.trav_id();
  };

  void mark_todo( node_index_t const& n )
  {
    ntk_.set_visited( n, ntk_.trav_id() );
  };

  /*! \brief Computes the initial required times.
   *
   * Traverse the network in reverse topological order and computes the reequired']
   * time of each node from its fanouts' required times.
   */
  void compute_required_times()
  {
    times_.reset( std::numeric_limits<double>::max() );

    uint32_t i = 0;
    ntk_.foreach_po( [&]( auto f ) {
      times_[f] = output_[i++];
    } );

    topo_sort_.foreach_node_reverse( [&]( node_index_t const& n ) {
      ntk_.foreach_output( n, [&]( auto const& fo ) {
        ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
          auto ni = ntk_.get_node( fi );
          update_required_time( fo, ii );
        } );
      } );
    } );
  }

  /*! \brief Update the required time of a signal's fanin from its required time.
   *
   * \param f A signal whose required time is known
   * \param i The index of it's fanin whose required time we want to update.
   */
  void update_required_time( signal_t const& f, uint32_t i )
  {
    node_index_t const n = ntk_.get_node( f );
    auto const& children = ntk_.get_children( n );

    double new_time = times_[f] - ntk_.get_max_pin_delay( f, i );
    if ( new_time < times_[children[i]] )
    {
      times_[children[i]] = new_time;
    }
  }

  /*! \brief Check if the required time at a node differs from the stored value
   */
  bool has_required_time_update( node_index_t n )
  {
    bool update = false;
    ntk_.foreach_output( n, [&]( auto const& f ) {
      update |= has_required_time_update( f );
    } );
    return update;
  }

  /*! \brief Check if the required time at a node differs from the stored value
   */
  bool has_required_time_update( signal_t const& f )
  {
    if ( ntk_.is_po( f ) )
    {
      return times_[f] != output_[ntk_.po_index( f )];
    }
    bool update = false;
    double new_time = std::numeric_limits<double>::max();
    ntk_.foreach_fanout( f, [&]( auto const& no ) {
      ntk_.foreach_output( no, [&]( auto const& fo ) {
        ntk_.foreach_fanin( no, [&]( auto const& fi, auto const ii ) {
          if ( fi == f )
          {
            new_time = std::min( times_[fo] - ntk_.get_max_pin_delay( fo, ii ), new_time );
          }
        } );
      } );
    } );
    return std::abs( new_time - times_[f] ) > std::numeric_limits<double>::epsilon();
  }

  void update_required_time( node_index_t n )
  {
    ntk_.foreach_output( n, [&]( auto const& f ) {
      update_required_time( f );
    } );
  }

  void update_required_time( signal_t const& f )
  {
    if ( ntk_.is_po( f ) )
    {
      times_[f] = output_[ntk_.po_index( f )];
      return;
    }
    double new_time = std::numeric_limits<double>::max();
    ntk_.foreach_fanout( f, [&]( auto const& no ) {
      ntk_.foreach_output( no, [&]( auto const& fo ) {
        ntk_.foreach_fanin( no, [&]( auto const& fi, auto const ii ) {
          if ( fi == f )
          {
            new_time = std::min( times_[fo] - ntk_.get_max_pin_delay( fo, ii ), new_time );
          }
        } );
      } );
    } );
    times_[f] = new_time;
  }

  void reset_required_tfi( node_index_t const& n )
  {
    auto f = ntk_.make_signal( n );
    return reset_required_tfi( f );
  }

  void reset_required_tfi( signal_t const& f )
  {
    auto const n = ntk_.get_node( f );
    if ( is_marked_todo( n ) )
      return;

    mark_todo( n );

    if ( !ntk_.is_po( f ) )
      times_[f] = infinite_time;
    else
      times_[f] = output_[ntk_.po_index( f )];

    if ( ntk_.is_pi( ntk_.get_node( f ) ) )
    {
      return;
    }

    ntk_.foreach_fanin( ntk_.get_node( f ), [&]( auto const& fi, auto ii ) {
      (void)ii;
      reset_required_tfi( fi );
    } );
  }

  void update_required( uint32_t last_level )
  {
    topo_sort_.foreach_node_reverse( last_level, [&]( auto const& n ) {
      if ( is_marked_todo( n ) )
        update_required_time( n );
    } );
  }

#pragma endregion

private:
  Ntk& ntk_;
  /* maintains the topological order of the network up to date */
  topo_sort_tracker<Ntk> topo_sort_;
  network::incomplete_signal_map<double, Ntk> times_;
  /* required times at the outputs */
  std::vector<double> output_;
  /* events */
  std::shared_ptr<typename mockturtle::network_events<Ntk>::add_event_type> add_event_;
  std::shared_ptr<typename mockturtle::network_events<Ntk>::modified_event_type> modified_event_;
};
} // namespace trackers

} // namespace analyzers

} // namespace rinox