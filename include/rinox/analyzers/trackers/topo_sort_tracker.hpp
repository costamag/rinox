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
 * \file topo_sort_tracker.hpp
 * \brief Extracts the topological order of a network and updates it
 *
 * \author Andrea Costamagna
 */

#pragma once

#include "../../network/tfo_manager.hpp"
#include <limits>
#include <mockturtle/utils/node_map.hpp>

namespace rinox
{

namespace analyzers
{

namespace trackers
{

/*! \brief Engine to efficiently maintain a topological order of a network
 *
 * This engine stores the nodes in topological order. It is based on the simple
 * observation that two nodes having the same shortest path to the PIs ( in
 * terms of the number of nodes ) cannot be one in the TFO of the other. Hence,
 * by classifying each node based on the shortest path to the PIs, it is easy
 * to maintain a topological order of the network at any time, with the same
 * algorithmic structure of a depth tracker.
 *
 * Nodes at the same depth are stored in a linked list, providing efficient
 * insertions/removals and making it easier to iterate in order.
 *
 * It listens to add, delete, and modify events from the network and updates its
 * internal topological structure incrementally rather than recomputing from scratch.
 *
 * Nodes in the same depth class are connected in the form of a linked list.
 *
 * \tparam Ntk the network type to be analyzed.
 *
 * \verbatim embed:rst

   Example

   .. code-block:: c++

      bound_network ntk( gates );
      auto const a = ntk.create_pi();
      auto const f1 = ntk.create_node( { a }, 0 );
      topo_sort_tracker topo_sort( ntk );
      auto const f2 = ntk.create_node( { f1 }, 0 );
   \endverbatim
 */
template<class Ntk>
class topo_sort_tracker
{
public:
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;
  node_index_t null = std::numeric_limits<node_index_t>::max();

private:
  struct node_info_t
  {
    node_index_t prev;
    node_index_t next;
    uint32_t level;
  };

public:
  topo_sort_tracker( Ntk& ntk )
      : ntk_( ntk ),
        nodes_( ntk ),
        tfo_( ntk )
  {
    init();
  }

  void init()
  {
    /* check if the network implements the needed functionalities */
    static_assert( mockturtle::is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( mockturtle::has_incr_trav_id_v<Ntk>, "Ntk does not implement the incr_trav_id method" );
    static_assert( mockturtle::has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
    static_assert( mockturtle::has_foreach_fanout_v<Ntk>, "Ntk does not implement the foreach_fanout method" );
    static_assert( mockturtle::has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
    static_assert( mockturtle::has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
    static_assert( mockturtle::has_set_value_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( mockturtle::has_value_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( mockturtle::has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( mockturtle::has_visited_v<Ntk>, "Ntk does not implement the visited method" );
    static_assert( mockturtle::has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
    static_assert( mockturtle::has_size_v<Ntk>, "Ntk does not implement the size method" );

    tails_.reserve( 200 );
    heads_.reserve( 200 );
    nodes_.resize();

    compute_topo_sort();

    add_event_ = ntk_.events().register_add_event( [&]( const node_index_t& n ) {
      nodes_.resize();
      auto level = compute_level( n );
      if ( level < heads_.size() )
      {
        node_index_t old_head = heads_[level];
        heads_[level] = n;
        nodes_[n] = { old_head, null, level };
        nodes_[old_head].next = n;
      }
      else
      {
        assert( level == heads_.size() );
        heads_.push_back( n );
        tails_.push_back( n );
        nodes_[n] = { null, null, level };
      }
    } );

    delete_event_ = ntk_.events().register_delete_event( [&]( const node_index_t& n ) {
      auto info = nodes_[n];
      node_index_t const prev = info.prev;
      node_index_t const next = info.next;
      uint32_t const level = info.level;
      if ( !is_head( n ) && !is_tail( n ) )
      {
        nodes_[prev].next = next;
        nodes_[next].prev = prev;
      }
      if ( is_head( n ) )
      {
        heads_[level] = nodes_[n].prev;
        if ( heads_[level] != null )
          nodes_[heads_[level]].next = null;
      }
      else if ( is_tail( n ) )
      {
        tails_[level] = nodes_[n].next;
        if ( tails_[level] != null )
          nodes_[tails_[level]].prev = null;
      }
    } );

    modified_event_ = ntk_.events().register_modified_event( [&]( const auto& n, auto old_children ) {
      /* The tfo of the new node is modified as it acquires new fanouts */
      update_topo_sort_tfo( n );

      /* The tfo of the old children is modified since they lose the old node */
      for ( auto const& f : old_children )
      {
        update_topo_sort_tfo( ntk_.get_node( f ) );
      }
    } );
  }

  ~topo_sort_tracker()
  {
    if ( add_event_ )
    {
      ntk_.events().release_add_event( add_event_ );
    }

    if ( delete_event_ )
    {
      ntk_.events().release_delete_event( delete_event_ );
    }

    if ( modified_event_ )
    {
      ntk_.events().release_modified_event( modified_event_ );
    }
  }

#pragma region Iterators
  template<typename Fn>
  void foreach_node( Fn&& fn )
  {
    for ( auto i = 0u; i < tails_.size(); ++i )
    {
      node_index_t tail = tails_[i];
      node_index_t head = tail;
      while ( head != null )
      {
        fn( head );
        head = nodes_[head].next;
      }
    }
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn )
  {
    for ( auto i = 0u; i < tails_.size(); ++i )
    {
      node_index_t tail = tails_[i];
      node_index_t head = tail;
      while ( head != null )
      {
        if ( !ntk_.is_pi( head ) && !ntk_.is_constant( head ) )
          fn( head );
        head = nodes_[head].next;
      }
    }
  }

  template<typename Fn>
  void foreach_node_reverse( Fn&& fn )
  {
    node_index_t head;

    for ( int i = heads_.size() - 1; i >= 0; --i )
    {
      node_index_t head = heads_[i];
      node_index_t tail = head;
      while ( tail != null )
      {
        fn( tail );
        tail = nodes_[tail].prev;
      }
    }
  }

  template<typename Fn>
  void foreach_node_reverse( uint32_t last_level, Fn&& fn )
  {
    node_index_t head;

    for ( int i = last_level; i >= 0; --i )
    {
      node_index_t head = heads_[i];
      node_index_t tail = head;
      while ( tail != null )
      {
        fn( tail );
        tail = nodes_[tail].prev;
      }
    }
  }
#pragma endregion

#pragma region Getters
public:
  [[nodiscard]] uint32_t get_level( node_index_t const& n ) const
  {
    return nodes_[n].level;
  }

  [[nodiscard]] std::vector<node_index_t> get_topological_order()
  {
    std::vector<node_index_t> order;
    foreach_node( [&]( auto const& n ) {
      order.push_back( n );
    } );
    return order;
  }

  [[nodiscard]] std::vector<node_index_t> get_reverse_order()
  {
    std::vector<node_index_t> order;
    foreach_node_reverse( [&]( auto const& n ) {
      order.push_back( n );
    } );
    return order;
  }
#pragma endregion

#pragma region Implementation details
private:
  bool is_head( node_index_t n )
  {
    auto const info = nodes_[n];
    return heads_[info.level] == n;
  }

  bool is_tail( node_index_t n )
  {
    auto const info = nodes_[n];
    return tails_[info.level] == n;
  }

  bool is_marked_ready( node_index_t const& n )
  {
    return ntk_.value( n ) == ntk_.trav_id();
  };

  void make_ready( node_index_t const& n )
  {
    ntk_.set_value( n, ntk_.trav_id() );
  };

  void compute_topo_sort()
  {
    nodes_.reset();
    if ( ntk_.num_pis() == 0 )
      return;

    ntk_.incr_trav_id();
    /* store the level 0 tail */
    node_index_t n = ntk_.pi_at( 0 );
    tails_ = { n };
    heads_ = { n };
    nodes_[n] = { null, null, 0 };
    for ( auto i = 1u; i < ntk_.num_pis(); ++i )
    {
      n = ntk_.pi_at( i );
      nodes_[n] = { heads_[0], null, 0 };
      nodes_[heads_[0]].next = n;
      heads_[0] = n;
      make_ready( n );
    }

    ntk_.foreach_po( [&]( auto const& f ) {
      compute_topo_sort_tfi( f );
    } );
  }

  /*! \brief Compute the depth times of the nodes in the TFI of a signal's node*/
  void compute_topo_sort_tfi( signal_t const& f )
  {
    node_index_t n = ntk_.get_node( f );
    if ( is_marked_ready( n ) || ntk_.is_pi( n ) )
    {
      return;
    }
    make_ready( n );

    ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
      compute_topo_sort_tfi( fi );
    } );
    auto level = compute_level( n );
    if ( level < tails_.size() )
    {
      nodes_[n] = { heads_[level], null, level };
      nodes_[heads_[level]].next = n;
      heads_[level] = n;
    }
    else
    {
      assert( level == tails_.size() );
      tails_.push_back( n );
      heads_.push_back( n );
      nodes_[n] = { null, null, level };
    }
  }

  /*! \brief Efficient update of the depth times in the TFO of a node.\
   */
  void update_topo_sort_tfo( node_index_t const& n )
  {
    /* mark all the nodes that must be considered */
    tfo_.init( n );

    std::vector<node_index_t> old_worklist{ n };
    std::vector<node_index_t> new_worklist;
    old_worklist.reserve( 100 );
    new_worklist.reserve( 100 );
    tfo_.mark_seen( n );

    bool progress{ true };
    /* iterate over all the nodes having some inputs ready */
    while ( progress && !old_worklist.empty() )
    {
      progress = false;
      for ( node_index_t const& u : old_worklist )
      {
        bool ready = true;
        ntk_.foreach_fanin( u, [&]( auto const& f, auto i ) {
          node_index_t ni = ntk_.get_node( f );
          ready &= ntk_.is_pi( ni ) || !tfo_.belongs_to_tfo( ni ) || tfo_.is_marked_ready( ni );
        } );
        tfo_.mark_seen( u );
        if ( ready )
        {
          progress = true;
          tfo_.mark_ready( u );
          double old_level = nodes_[u].level;
          auto new_level = compute_level( u );
          if ( new_level != old_level )
          {
            ntk_.foreach_fanout( u, [&]( node_index_t const& o ) {
              if ( !tfo_.is_marked_seen( o ) ) // only insert if not seen
              {
                new_worklist.push_back( o );
                tfo_.mark_seen( o );
              }
            } );

            if ( !is_tail( u ) && !is_head( u ) )
            {
              nodes_[nodes_[u].prev].next = nodes_[u].next;
              nodes_[nodes_[u].next].prev = nodes_[u].prev;
            }
            else if ( is_tail( u )&&!is_head(u) ) // but not head
            {
              node_index_t new_tail = nodes_[u].next;
              tails_[old_level] = new_tail;
              nodes_[new_tail].prev = null;
            }
            else if ( !is_tail( u )&& is_head(u) )
            {
              node_index_t new_head = nodes_[u].prev;
              heads_[old_level] = new_head;
              nodes_[new_head].next = null;
            }
            else
            {
              tails_.erase( tails_.begin()+old_level);
              heads_.erase( heads_.begin()+old_level);
            }

            if ( new_level < heads_.size() )
            {
              nodes_[u] = { heads_[new_level], null, new_level };
              nodes_[heads_[new_level]].next = u;
              heads_[new_level] = u;
            }
            else
            {
              assert( new_level == heads_.size() );
              nodes_[u] = { null, null, new_level };
              heads_.push_back( u );
              tails_.push_back( u );
            }
          }
        }
        else
        {
          new_worklist.push_back( u );
        }
      }
      old_worklist.swap( new_worklist );
      new_worklist.clear();
    }
//    if ( !progress )
//    {
//      std::cerr << "[e] Infinite loop in level times analyzer" << std::endl;
//    }
  }

  [[nodiscard]] uint32_t compute_level( node_index_t const& n )
  {
    if ( ntk_.is_pi( n ) )
    {
      return 0;
    }
    else
    {
      auto level = std::numeric_limits<uint32_t>::min();
      ntk_.foreach_fanin( n, [&]( auto const& fi, auto const ii ) {
        node_index_t ni = ntk_.get_node( fi );
        level = std::max( level, nodes_[ni].level + 1 );
      } );
      return level;
    }
  }
#pragma endregion

private:
  Ntk& ntk_;
  /* container of the shortest path to a node and its links to nodes on the same level */
  mockturtle::incomplete_node_map<node_info_t, Ntk> nodes_;
  /* first node in each depth class */
  std::vector<node_index_t> tails_;
  std::vector<node_index_t> heads_;
  network::tfo_manager<Ntk> tfo_;
  /* events */
  std::shared_ptr<typename mockturtle::network_events<Ntk>::add_event_type> add_event_;
  std::shared_ptr<typename mockturtle::network_events<Ntk>::delete_event_type> delete_event_;
  std::shared_ptr<typename mockturtle::network_events<Ntk>::modified_event_type> modified_event_;
};
} // namespace trackers

} // namespace analyzers

} // namespace rinox