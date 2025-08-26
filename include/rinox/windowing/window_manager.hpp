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
  \file window_manager.hpp
  \brief Construction of windows for mapped networks.

  \author Andrea Costamagna
*/

#pragma once

#include <mockturtle/utils/node_map.hpp>
#include "../network/visualize.hpp"

namespace rinox
{

namespace windowing
{

template<class Ntk>
struct window_t
{
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;

  node_index_t pivot;
  std::vector<node_index_t> tfos;
  std::vector<node_index_t> mffc;
  std::vector<signal_t> divs;
  std::vector<signal_t> outputs;
  std::vector<signal_t> inputs;
};

struct default_window_manager_params
{
  static constexpr uint32_t max_num_leaves = 8;
  uint32_t max_num_divisors = 128;
  bool preserve_depth = true;
  int32_t odc_levels = 0;
  uint32_t skip_fanout_limit_for_divisors{ 100 };
};

struct window_manager_stats
{
  bool valid = false;
};

template<class Ntk, typename Params = default_window_manager_params>
class window_manager
{
public:
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;

public:
  window_manager( Ntk& ntk, Params const& ps, window_manager_stats& st )
      : ntk_( ntk ),
        color_map_( ntk ),
        ps_( ps ),
        st_( st )
  {
  }
#if 0
  bool run( node_index_t const& n )
  {
    st_.valid = false;
    init( n );

    // label the MFFC nodes
    collect_mffc_nodes();
    for ( auto const& m : window_.mffc )
    {
      ntk_.set_visited( m, ntk_.trav_id() );
      make_alien( m );
    }
    window_.mffc.clear();

    window_.inputs.clear();
    ntk_.foreach_output( n, [&]( auto const& f ) {
      window_.inputs.push_back( f );
      window_.divs.push_back( f );
    } );

    // expand the leaves until reaching the boundary of the MFFC
    expand_leaves( [&]( node_index_t const& v ) { return ( ntk_.visited( v ) == ntk_.trav_id() ) && !ntk_.is_pi( v ); },
                   [&]( node_index_t const& v ) { make_mffc( v ); } );

    collect_mffc_nodes();

    window_.divs = window_.inputs;
    /* expand toward the tfo */
    if ( ps_.odc_levels > 0 )
    {
      collect_nodes_tfo();
      // expand the divisors set to try removing upper-leaves with reconvergence
      collect_side_divisors();
    }
    else
    {
      ntk_.foreach_output( n, [&]( auto f ) {
        window_.outputs.push_back( f );
      } );
    }

    // expand the leaves to find reconvergences
    expand_leaves( [&]( node_index_t const& n ) { return !ntk_.is_pi( n ); },
                   [&]( node_index_t const& v ) { make_divisor( v ); } );

    // expand the divisors set to try removing upper-leaves with reconvergence
    collect_side_divisors();

    if ( ps_.odc_levels > 0 )
    {
      topological_sort( window_.tfos );
      topological_sort( window_.outputs );
    }
    topological_sort( window_.mffc );
    topological_sort( window_.divs );
    topological_sort( window_.inputs );

    st_.valid = true;
    return true;
  }
#endif

  bool run( node_index_t const& n )
  {
    st_.valid = false;
    init( n );

    ntk_.incr_trav_id();

    window_.pivot = n;

    if ( ntk_.fanin_size( n ) == 0 )
      return false;
    // after this the mffc nodes are collected and marked
    // and the first nodes not in the mffc are considered as leaves
    // marked as such, and collected under inputs
    collect_mffc_nodes();

    /* expand toward the tfo
       mark the tfo nodes, collect the output signals, and add the fanins not yet included as inputs*/
    collect_tfos_nodes();

    // expand the leaves to find reconvergences
    bool expanded = true;
    while ( expanded )
    {
      expanded = false;
      expanded &= collect_leaf_nodes();
      expanded &= collect_divs_nodes();
    }

    topological_sort( window_.tfos );
    topological_sort( window_.outputs );
    topological_sort( window_.mffc );
    topological_sort( window_.divs );
    topological_sort( window_.inputs );

    st_.valid = true;
    st_.valid &= window_.divs.size() <= ps_.max_num_divisors;
    st_.valid &= window_.inputs.size() <= ps_.max_num_leaves;
    return st_.valid;
  }

  std::vector<signal_t> const& get_divisors() const
  {
    return window_.divs;
  }

  signal_t const& get_divisor( uint32_t const& index ) const
  {
    return window_.divs[index];
  }

  std::vector<node_index_t> const& get_tfos() const
  {
    return window_.tfos;
  }

  std::vector<signal_t> const& get_outputs() const
  {
    return window_.outputs;
  }

  std::vector<signal_t> const& get_inputs() const
  {
    return window_.inputs;
  }

  std::vector<signal_t> const& get_divs() const
  {
    return window_.divs;
  }

  std::vector<node_index_t> const& get_mffc() const
  {
    return window_.mffc;
  }

  node_index_t const& get_pivot() const
  {
    return window_.pivot;
  }

  window_t<Ntk> get_window() const
  {
    return window_;
  }

#pragma region Miscellanea
private:
  void init( node_index_t const& n )
  {
    color_map_.resize();
    color_++;
    ntk_.incr_trav_id();
    window_.outputs.clear();
    window_.tfos.clear();
    window_.mffc.clear();
    window_.divs.clear();
    window_.inputs.clear();
  }

  void topological_sort( std::vector<node_index_t>& nodes )
  {
    std::stable_sort( nodes.begin(), nodes.end(), [&]( auto const& a, auto const& b ) {
      return ntk_.level( a ) < ntk_.level( b );
    } );
  }

  void topological_sort( std::vector<signal_t>& signals )
  {
    std::stable_sort( signals.begin(), signals.end(), [&]( auto const& a, auto const& b ) {
      return ntk_.level( ntk_.get_node( a ) ) < ntk_.level( ntk_.get_node( b ) );
    } );
  }

public:
  bool is_tfo( node_index_t const& n ) const
  {
    return is_contained( n ) && ( color_map_[n] == ( 4u | ( color_ << 3u ) ) );
  }

  void make_tfo( node_index_t const& n )
  {
    color_map_[n] = ( 4u | ( color_ << 3u ) );
  }

  bool is_contained( node_index_t const& n ) const
  {
    return (n < color_map_.size()) && color_map_.has( n ) && ( color_map_[n] >> 3u == color_ );
  }

  void mark_contained()
  {
    ntk_.incr_trav_id();
    for ( auto const& f : window_.divs )
    {
      auto const n = ntk_.get_node( f );
      make_divisor( n );
    }

    for ( auto const& f : window_.inputs )
    {
      auto const n = ntk_.get_node( f );
      make_input( n );
    }

    for ( auto const& n : window_.tfos )
    {
      make_tfo( n );
    }

    for ( auto const& f : window_.outputs )
    {
      auto const n = ntk_.get_node( f );
      make_output( n );
    }

    for ( auto const& n : window_.mffc )
    {
      make_mffc( n );
    }
  }

  bool is_output( node_index_t const& n ) const
  {
    return is_contained( n ) && ( color_map_[n] == ( 5u | ( color_ << 3u ) ) );
  }

  void make_output( node_index_t const& n )
  {
    color_map_[n] = ( 5u | ( color_ << 3u ) );
  }

  void make_alien( node_index_t const& n )
  {
    color_map_[n] = 0;
  }

  bool is_input( node_index_t const& n ) const
  {
    return is_contained( n ) && ( color_map_[n] == ( 3u | ( color_ << 3u ) ) );
  }

  void make_input( node_index_t const& n )
  {
    color_map_[n] = ( 3u | ( color_ << 3u ) );
  }

  bool is_mffc( node_index_t const& n ) const
  {
    return is_contained( n ) && ( color_map_[n] == ( 1u | ( color_ << 3u ) ) );
  }

  void make_mffc( node_index_t const& n )
  {
    color_map_[n] = ( 1u | ( color_ << 3u ) );
  }

  bool is_divisor( node_index_t const& n ) const
  {
    return is_contained( n ) && ( color_map_[n] == ( 2u | ( color_ << 3u ) ) );
  }

  void make_divisor( node_index_t const& n )
  {
    color_map_[n] = ( 2u | ( color_ << 3u ) );
  }

  bool is_valid() const
  {
    return st_.valid;
  }
#pragma endregion

#pragma region TFO
  void collect_tfos_nodes()
  {
    ntk_.foreach_output( window_.pivot, [&]( auto f ) {
      if ( ntk_.is_po( f ) || ( ps_.odc_levels <= 0 ) )
        window_.outputs.push_back( f );
    } );
    if ( ps_.odc_levels <= 0 )
    {
      return;
    }

    std::vector<node_index_t> fronteer;
    ntk_.foreach_fanout( window_.pivot, [&]( auto const& no ) {
      if ( ntk_.visited( no ) != ntk_.trav_id() )
      {
        ntk_.set_visited( no, ntk_.trav_id() );
        fronteer.push_back( no );
      }
    } );

    uint32_t level = 0;
    while ( fronteer.size() > 0 && ( level < ps_.odc_levels ) )
    {
      std::vector<node_index_t> new_fronteer;
      level++;
      // add the POs to the outputs list
      for ( auto const& n : fronteer )
      {
        window_.tfos.push_back( n );
        make_tfo( n );
        ntk_.foreach_output( n, [&]( auto const& f ) {
          if ( ntk_.is_po( f ) )
          {
            window_.outputs.push_back( f );
          }
          if ( level == ps_.odc_levels )
          {
            bool dep = true;
            ntk_.foreach_fanout( f, [&]( auto const& no ) {
              dep &= is_contained( no );
            } );
            if ( !dep )
              window_.outputs.push_back( f );
          }
        } );
      }
      // expand the fronteer
      for ( auto const& n : fronteer )
      {
        ntk_.foreach_fanout( n, [&]( auto const& no ) {
          if ( !is_tfo( no ) )
          {
            make_tfo( no );
            new_fronteer.push_back( no );
          }
        } );
      }
      for ( auto const& n : fronteer )
      {
        ntk_.foreach_fanin( n, [&]( auto const& fi ) {
          auto const ni = ntk_.get_node( fi );
          if ( !is_contained( ni ) && ( ntk_.visited( ni ) != ntk_.trav_id() ) )
          {
            window_.inputs.push_back( fi );
            ntk_.set_visited( ni, ntk_.trav_id() );
            make_input( ni );
          }
        } );
      }
      fronteer = new_fronteer;
    }
    window_.inputs.erase( std::remove_if( window_.inputs.begin(), window_.inputs.end(), [&]( auto const& f ) {
                            auto const ni = ntk_.get_node( f );
                            return !is_input( ni );
                          } ),
                          window_.inputs.end() );
    window_.divs = window_.inputs;
  }
#pragma endregion

  bool collect_divs_nodes()
  {
    int cnt = 0;
    bool done = false, expanded = false;
    while ( !done )
    {
      done = true;
      auto const num_divs = window_.divs.size();
      for ( auto i = 0u; i < num_divs; ++i )
      {
        ntk_.foreach_fanout( window_.divs[i], [&]( auto const no ) {
          if ( window_.divs.size() >= ps_.max_num_divisors )
            return;

          if ( is_input( no ) )
          {
            bool all_contained = true;
            ntk_.foreach_fanin( no, [&]( auto const& fi, auto const ii ){
              all_contained &= is_contained( ntk_.get_node( fi ) );
            } );
            if ( all_contained )
              make_divisor( no );
          }
          else if ( !is_contained( no ) && ( ntk_.visited( no ) != ntk_.trav_id() ) )
          {
            bool in_divs = true;
            ntk_.foreach_fanin( no, [&]( auto const& fi ) {
              auto const ni = ntk_.get_node( fi );
              in_divs &= is_divisor( ni ) || is_input( ni );
            } );
            if ( in_divs && ( ( ntk_.num_outputs( no ) + window_.divs.size() ) < ps_.max_num_divisors ) )
            {
              ntk_.foreach_output( no, [&]( auto const& fo ) {
                window_.divs.push_back( fo );
                make_divisor( no );
                ntk_.set_visited( no, ntk_.trav_id() );
              } );
              done = false;
              expanded = true;
            }
          }
        } );
        if ( window_.divs.size() >= ps_.max_num_divisors )
          done = true;
      }

      window_.inputs.erase( std::remove_if( window_.inputs.begin(), window_.inputs.end(), [&]( auto const& f ) {
                              auto const ni = ntk_.get_node( f );
                              return !is_input( ni );
                            } ),
                            window_.inputs.end() );
    }
    return expanded;
  }

  bool collect_leaf_nodes()
  {
    bool expanded = false, done = false;
    while ( !done )
    {
      done = true;
      int best_cost = ps_.max_num_leaves - window_.inputs.size() + 1;
      std::optional<node_index_t> best_node = std::nullopt;
      auto const num_inputs = window_.inputs.size();
      for ( auto i = 0u; i < num_inputs; ++i )
      {
        auto const l = window_.inputs[i];
        node_index_t leaf = ntk_.get_node( l );
        int cost = compute_leaf_cost( leaf );
        if ( cost < best_cost )
        {
          best_cost = cost;
          best_node = std::make_optional( leaf );
        }
      }
      if ( best_node )
      {
        ntk_.foreach_fanin( *best_node, [&]( auto const& fi, auto ii ) {
          auto const ni = ntk_.get_node( fi );
          if ( !is_contained( ni ) )
          {
            ntk_.foreach_output( ni, [&]( auto const& fo ) {
              window_.inputs.push_back( fo );
              window_.divs.push_back( fo );
            } );
            make_input( ni );
            expanded = true;
            done = false;
          }
        } );
        make_divisor( *best_node );
      }
      window_.inputs.erase( std::remove_if( window_.inputs.begin(), window_.inputs.end(), [&]( auto const& f ) {
                              auto const ni = ntk_.get_node( f );
                              return !is_input( ni );
                            } ),
                            window_.inputs.end() );
    }
    return expanded;
  }

#pragma region Leaves
  int compute_leaf_cost( node_index_t const& n )
  {
    if ( ntk_.is_pi( n ) )
      return ps_.max_num_leaves;

    int cost = -static_cast<int>( ntk_.num_outputs( n ) );
    ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
      auto const ni = ntk_.get_node( fi );
      if ( !is_contained( ni ) )
        cost += ntk_.num_outputs( ni );
    } );
    return cost;
  }
#pragma endregion Leaves

#pragma region Mffc
  void collect_mffc_nodes()
  {
    make_mffc( window_.pivot );
    window_.mffc = { window_.pivot };
    ntk_.set_visited( window_.pivot, ntk_.trav_id() );

    std::vector<node_index_t> fronteer;
    ntk_.foreach_fanin( window_.pivot, [&]( auto const& fi ) {
      auto const ni = ntk_.get_node( fi );
      if ( ntk_.visited( ni ) != ntk_.trav_id() )
      {
        ntk_.set_visited( ni, ntk_.trav_id() );
        fronteer.push_back( ni );
      }
    } );

    while ( fronteer.size() > 0 )
    {
      std::vector<node_index_t> new_fronteer;
      for ( auto const& n : fronteer )
      {
        bool in_mffc = !ntk_.is_pi( n ) && !is_contained( n ) && !ntk_.is_po( n );
        if ( in_mffc )
        {
          ntk_.foreach_fanout( n, [&]( auto const& no ) {
            in_mffc &= is_mffc( no );
          } );
        }
        if ( in_mffc )
        {
          make_mffc( n );
          window_.mffc.push_back( n );
          ntk_.foreach_fanin( n, [&]( auto const& fi ) {
            auto const ni = ntk_.get_node( fi );
            if ( ntk_.visited( ni ) != ntk_.trav_id() )
            {
              new_fronteer.push_back( ni );
              ntk_.set_visited( ni, ntk_.trav_id() );
            }
          } );
        }
        else if ( !is_contained( n ) )
        {
          ntk_.foreach_output( n, [&]( auto const& f ) {
            window_.inputs.push_back( f );
            window_.divs.push_back( f );
          } );
          make_input( n );
        }
      }
      fronteer = new_fronteer;
    }
  }

  /* ! \brief Dereference the node's MFFC */
  void node_deref_rec( node_index_t const& n )
  {
    if ( ntk_.is_pi( n ) )
      return;
    ntk_.foreach_fanin( n, [&]( const auto& f ) {
      auto const& p = ntk_.get_node( f );
      if ( ntk_.is_pi( p ) )
        return;
      ntk_.decr_fanout_size( p );
      if ( ntk_.fanout_size( p ) == 0 )
      {
        make_mffc( p );
        window_.mffc.push_back( p );
        node_deref_rec( p );
      }
    } );
  }

  /* ! \brief Reference the node's MFFC */
  void node_ref_rec( node_index_t const& n )
  {
    if ( ntk_.is_pi( n ) )
      return;
    ntk_.foreach_fanin( n, [&]( const auto& f ) {
      auto const& p = ntk_.get_node( f );
      if ( ntk_.is_pi( p ) )
        return;
      auto v = ntk_.fanout_size( p );
      ntk_.incr_fanout_size( p );
      if ( v == 0 && !ntk_.is_pi( p ) )
      {
        node_ref_rec( p );
      }
    } );
  }
#pragma endregion

#pragma region Divisors
  void collect_side_divisors()
  {
    uint32_t max_level = std::numeric_limits<uint32_t>::min();
    for ( auto const& f : window_.outputs )
      max_level = std::max( max_level, ntk_.level( ntk_.get_node( f ) ) );

    bool done = false;
    while ( !done )
    {
      done = true;
      // check if any leaf is contained by other leaves
      for ( auto const& f : window_.inputs )
      {
        auto const n = ntk_.get_node( f );
        if ( is_input( n ) && !ntk_.is_pi( n ) )
        {
          bool contained = true;
          ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
            auto const ni = ntk_.get_node( fi );
            contained &= is_contained( ni ) && !is_tfo( ni ) && !is_mffc( ni ) && !is_output( ni );
          } );
          if ( contained )
          {
            make_divisor( n );
          }
        }
      }

      auto& inputs = window_.inputs;
      inputs.erase( std::remove_if( inputs.begin(),
                                    inputs.end(),
                                    [&]( signal_t const& f ) { return is_divisor( ntk_.get_node( f ) ); } ),
                    inputs.end() );

      std::vector<signal_t> new_divs;
      for ( auto const& f : window_.divs )
      {
        auto const n = ntk_.get_node( f );
        ntk_.foreach_fanout( n, [&]( auto const& no ) {
          if ( is_contained( no ) )
            return;
          if ( ps_.preserve_depth && ( ntk_.level( no ) >= max_level ) )
            return;
          bool add = true;
          ntk_.foreach_fanin( no, [&]( auto const& fi, auto ii ) {
            auto ni = ntk_.get_node( fi );
            add &= is_contained( ni ) && !is_tfo( ni ) && !is_mffc( ni ) && !is_output( ni );
          } );
          if ( add )
          {
            make_divisor( no );
            ntk_.foreach_output( no, [&]( auto const& fo ) {
              new_divs.push_back( fo );
            } );
            done = false;
          }
          return;
        } );
      }
      for ( auto d : new_divs )
        window_.divs.push_back( d );
    }
  }

public:
  size_t num_inputs() const
  {
    return window_.inputs.size();
  }

  size_t num_outputs() const
  {
    return window_.outputs.size();
  }

  size_t num_divisors() const
  {
    return window_.divs.size();
  }

  size_t size() const
  {
    return window_.divs.size() + window_.outputs.size() + ( window_.tfos.size() * Ntk::max_num_outputs );
  }

  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_input( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < window_.inputs.size(); i++ )
    {
      fn( window_.inputs[i], i );
    }
  }

  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_divisor( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < window_.divs.size(); i++ )
    {
      fn( window_.divs[i], i );
    }
  }

  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_mffc( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < window_.mffc.size(); i++ )
    {
      fn( window_.mffc[i], i );
    }
  }

  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_tfo( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < window_.tfos.size(); i++ )
    {
      fn( window_.tfos[i], i );
    }
  }

  /*! \brief iterator to apply a lambda to all the nodes */
  template<typename Fn>
  void foreach_output( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < window_.outputs.size(); i++ )
    {
      fn( window_.outputs[i], i );
    }
  }

  void print()
  {
    int i = 0;
    for ( auto const& l : window_.inputs )
      std::cout << "n[" << ntk_.get_node( l ) << "] = ntk.create_pi();" << std::endl;

    i = 0;
    for ( auto const& d : window_.divs )
    {
      if ( d.output == 0 )
      {
        auto const ids = ntk_.get_binding_ids( ntk_.get_node( d ) );
        std::cout << "n[" << ntk_.get_node( d ) << "] = ntk.create_node( { ";
        ntk_.foreach_fanin( d, [&]( auto const& fi, auto ii ){
          if ( ii < ( ntk_.fanin_size( d.index ) - 1 ) )
            std::cout << "signal{ n[" << fi.index << "], " << fi.output << " } , ";
          else
            std::cout << "signal{ n[" << fi.index << "], " << fi.output << " } } , { ";
        } );
        for ( auto i = 0u; i < (ids.size() - 1); ++i )
        {
          std::cout << ids[i] << ", ";
        }
        std::cout << ids.back() << " } ).index;" << std::endl;
      }
      std::cout << "ntk.create_po( { n[" << d.index << "], " << d.output << " } );" << std::endl;
    }

    for ( auto const& m : window_.mffc )
    {
      auto const ids = ntk_.get_binding_ids( m );
      std::cout << "n[" << m << "] = ntk.create_node( { ";
      ntk_.foreach_fanin( m, [&]( auto const& fi, auto ii ){
        if ( ii < ( ntk_.fanin_size( m ) - 1 ) )
          std::cout << "signal{ n[" << fi.index << "], " << fi.output << " } , ";
        else
          std::cout << "signal{ n[" << fi.index << "], " << fi.output << " } } , { ";
      } );
      for ( auto i = 0u; i < (ids.size() - 1); ++i )
      {
        std::cout << ids[i] << ", ";
      }
      std::cout << ids.back() << " } ).index;" << std::endl;
    }

    for ( auto const& m : window_.tfos )
    {
      auto const ids = ntk_.get_binding_ids( m );
      std::cout << "n[" << m << "] = ntk.create_node( ";
      ntk_.foreach_fanin( m, [&]( auto const& fi, auto ii ){
        if ( ii < ( ntk_.fanin_size( m ) - 1 ) )
          std::cout << "signal{ n[" << fi.index << "], " << fi.output << " } , ";
        else
          std::cout << "signal{ n[" << fi.index << "], " << fi.output << " } } , { ";
      } );
      for ( auto i = 0u; i < ids.size() - 1; ++i )
      {
        std::cout << ids[i] << ", ";
      }
      std::cout << ids.back() << " } ).index;" << std::endl;
    }

    i = 0;
    for ( auto const& l : window_.outputs )
      std::cout << "ntk.create_po( signal{ n[" << l.index << "], " << l.output << " } ); " << std::endl;
  }

private:
  Ntk& ntk_;
  window_t<Ntk> window_;
  mockturtle::incomplete_node_map<uint32_t, Ntk> color_map_;
  uint32_t color_ = 1u;
  Params const& ps_;
  window_manager_stats& st_;
};

} // namespace windowing

} // namespace rinox