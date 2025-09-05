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
 * \file tfo_manager.hpp
 * \brief Collects the TFO of a node and maintains flags on their exploration
 *
 * \author Andrea Costamagna
 */

#pragma once

#include <limits>
#include <mockturtle/utils/node_map.hpp>

namespace rinox
{

namespace network
{

/*! \brief Manager for the transitive fanout
 *
 * Data structure to extract and manipulate the TFO of a network.
 */
template<typename Ntk>
class tfo_manager
{
public:
  using node_index_t = typename Ntk::node;
  using index_t = uint64_t;
  using signal_t = typename Ntk::signal;
  /*! \brief Node information
   *
   * Contains flags to label a node based on:
   * - The belonging to the TFO of another node
   * - The readiness of its information ( arrival time )
   * - Whether it was seen or not in the process
   */
  struct node_info_t
  {
    node_info_t( uint64_t index, uint64_t ready, uint64_t seen )
        : index( index ), ready( ready ), seen( seen )
    {}

    node_info_t( uint64_t index )
        : node_info_t( index, 0u, 0u )
    {}

    node_info_t()
        : node_info_t( 0 ){}

      uint64_t index;
      uint64_t ready;
      uint64_t seen ;
  };

public:
  tfo_manager( Ntk& ntk )
      : ntk_( ntk ),
        map_( ntk )
  {
    // Force the map to allocate full network size with valid default values
  }

  /*! \brief Mark the nodes in the root's TFO */
  void init( node_index_t const& root )
  {
    map_.resize();
    root_ = root;
    assert( root_ < ntk_.size() );
    mark_tfo( root );
  }

  bool belongs_to_tfo( node_index_t const& n )
  {
    ensure_safe_access( n, "belongs_to_tfo" );
    auto const nindex = map_[n].index;
    return nindex == root_;
  }

  bool is_marked_ready( node_index_t const& n )
  {
    return map_[n].ready > 0u;
  }

  void mark_ready( node_index_t const& n )
  {
    map_[n].index = root_;
    map_[n].ready = 1;
  }

  bool is_marked_seen( node_index_t const& n )
  {
    return ntk_.is_pi( n ) || ( map_[n].seen > 0u );
  }

  void mark_seen( node_index_t const& n )
  {
    map_[n].index = root_;
    map_[n].seen = 1;
    map_[n].ready = 0;
  }

private:
  void mark_tfo( node_index_t const& n )
  {
    ensure_safe_access( n, "belongs_to_tfo" );

    if ( ntk_.is_pi( n ) || belongs_to_tfo( n ) || ntk_.is_dead( n ) )
      return;

    if ( !map_.has( n ) || !belongs_to_tfo( n ) )
    {
      make_tfo( n );
    }

    if ( ntk_.fanout_size( n ) == 0 )
      return;
    ntk_.foreach_fanout( n, [&]( auto u ) {
      mark_tfo( u );
    } );
  }

  void make_tfo( node_index_t const& n )
  {
    ensure_safe_access( n, "make_tfo" );
    map_[n] = node_info_t( root_, 0, 0 );
  }

  void ensure_map_capacity( node_index_t const& n )
  {
    auto idx = ntk_.node_to_index( n );
    while ( idx >= map_.size() )
    {
      map_.resize();
    }
  }

  void ensure_safe_access( node_index_t const& n, const char* caller )
  {
    auto idx = n;
    if ( idx >= map_.size() )
    {
      std::cerr << "[DEBUG] Out-of-bounds access in " << caller
                << ": node index = " << idx
                << ", map size = " << map_.size()
                << ", ntk_.size() = " << ntk_.size()
                << std::endl;
      std::abort(); // force a clean failure
    }
  }

private:
  /*! \brief Root node defining the TFO */
  node_index_t root_;
  /*! \brief Network where the TFO is analyzed */
  Ntk& ntk_;
  /*! \brief Container of the information of each node */
  mockturtle::incomplete_node_map<node_info_t, Ntk> map_;
};

} // namespace network

} // namespace rinox