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

#include "../math/math.hpp"
#include "../network/visualize.hpp"
#include <mockturtle/utils/node_map.hpp>

namespace rinox
{

namespace windowing
{

struct default_window_params_t
{
  static constexpr uint32_t tfos_capacity = 16u;
  static constexpr uint32_t mffc_capacity = 64u;
  static constexpr uint32_t divs_capacity = 1024u;
  static constexpr uint32_t outs_capacity = 16u;
  static constexpr uint32_t inps_capacity = 32u;
};

template<class Ntk, typename Params = default_window_params_t>
struct window_t
{
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;
  using static_params = Params;

  void clear()
  {
    tfos.clear();
    mffc.clear();
    divs.clear();
    outputs.clear();
    inputs.clear();
  }

  node_index_t pivot;
  rinox::math::static_vector<node_index_t, static_params::tfos_capacity> tfos;
  rinox::math::static_vector<node_index_t, static_params::mffc_capacity> mffc;
  rinox::math::static_vector<signal_t, static_params::divs_capacity> divs;
  rinox::math::static_vector<signal_t, static_params::outs_capacity> outputs;
  rinox::math::static_vector<signal_t, static_params::inps_capacity> inputs;
};

struct default_window_manager_params
{
  static constexpr uint32_t max_num_leaves = 8;
  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 10;
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
  static constexpr uint32_t tfos_capacity = window_t<Ntk>::static_params::tfos_capacity;
  static constexpr uint32_t outs_capacity = window_t<Ntk>::static_params::outs_capacity;
  static constexpr uint32_t mffc_capacity = window_t<Ntk>::static_params::mffc_capacity;
  static constexpr uint32_t divs_capacity = window_t<Ntk>::static_params::divs_capacity;
  static constexpr uint32_t inps_capacity = window_t<Ntk>::static_params::inps_capacity;
  using input_signals_t = rinox::math::static_vector<signal_t, inps_capacity>;
  using front_signals_t = rinox::math::static_vector<signal_t, inps_capacity + outs_capacity>;
  using front_nodes_t = rinox::math::static_vector<node_index_t, divs_capacity>;

public:
  window_manager( Ntk& ntk, Params const& ps, window_manager_stats& st )
      : ntk_( ntk ),
        color_map_( ntk_ ),
        ps_( ps ),
        st_( st )
  {
  }

  bool run( node_index_t const& n )
  {
    st_.valid = false;
    clear();

    if ( ntk_.fanin_size( n ) == 0 )
      return false;
    if ( ntk_.is_pi( n ) )
      return false;

    init( n );
    bool expanded = true;
    tfi_expansion();
    side_expansion();
    if ( ps_.odc_levels )
    {
      tfo_expansion();
      tfi_expansion();
      side_expansion();
    }

    rinox::graphics::visualize_multioutput_window( ntk_, window_, "win_" + std::to_string( n ) );
    remove_redundant();
    topological_sort<tfos_capacity>( window_.tfos );
    topological_sort<outs_capacity>( window_.outputs );
    topological_sort<mffc_capacity>( window_.mffc );
    topological_sort<divs_capacity>( window_.divs );
    topological_sort<inps_capacity>( window_.inputs );

    st_.valid = true;
    st_.valid &= window_.divs.size() <= ps_.max_num_divisors;
    st_.valid &= window_.inputs.size() <= ps_.max_num_leaves;

    return st_.valid;
  }

private:
  void remove_input( node_index_t const& n )
  {
    // erase signals corresponding to n from inputs
    for ( int i = window_.inputs.size() - 1; i >= 0; --i )
    {
      if ( ntk_.get_node( window_.inputs[i] ) == n )
        window_.inputs.erase( i );
    }

    // erase signals corresponding to n from divisors
    for ( int i = window_.divs.size() - 1; i >= 0; --i )
    {
      if ( ntk_.get_node( window_.divs[i] ) == n )
        window_.divs.erase( i );
    }

    // re-label node as alien before it becomes TFO
    make_alien( n );
  }

  void tfo_expansion()
  {
    front_nodes_.clear();
    front_.clear();

    // --- 0) Early budget check for pivot fanouts ---
    std::vector<node_index_t> pivot_fos;
    ntk_.foreach_fanout( window_.pivot, [&]( auto const& fo ) {
      pivot_fos.push_back( fo );
    } );

    std::set<node_index_t> new_inputs_all;
    for ( auto const& fo : pivot_fos )
    {
      ntk_.foreach_fanin( fo, [&]( auto const& fi, auto ) {
        auto const ni = ntk_.get_node( fi );
        if ( !is_contained( ni ) )
          new_inputs_all.insert( ni );
      } );
    }

    if ( window_.inputs.size() + new_inputs_all.size() > ps_.max_num_leaves )
    {
      // Pivot must remain sole output cut
      return;
    }

    // --- 1) Collect immediate fanouts as TFO ---
    for ( auto const& fo : pivot_fos )
    {
      make_tfo( fo );
      window_.tfos.insert_unique( fo );
      front_nodes_.push_back( fo );
    }

    // Add required inputs
    for ( auto ni : new_inputs_all )
    {
      if ( !is_contained( ni ) )
      {
        make_input( ni );
        ntk_.foreach_output( ni, [&]( auto const& sig ) {
          window_.inputs.insert_unique( sig );
          window_.divs.insert_unique( sig );
        } );
      }
    }

    // --- 2) Expand forward guided by cost ---
    while ( !front_nodes_.empty() )
    {
      auto tfo_opt = choose_removable_tfo();
      if ( !tfo_opt )
        break;
      auto const tfo = *tfo_opt;

      // Gather new fanins if we expand tfo
      std::set<node_index_t> new_inputs;
      ntk_.foreach_fanout( tfo, [&]( auto const& no ) {
        ntk_.foreach_fanin( no, [&]( auto const& fi, auto ) {
          auto const ni = ntk_.get_node( fi );
          if ( !is_contained( ni ) )
            new_inputs.insert( ni );
        } );
      } );

      // Budget check
      size_t budget_after = window_.inputs.size();
      if ( is_input( tfo ) )
        budget_after -= 1;
      for ( auto ni : new_inputs )
      {
        if ( !is_contained( ni ) )
          budget_after += 1;
      }

      if ( budget_after > ps_.max_num_leaves )
      {
        // Cannot expand → leave tfo as frontier
      }
      else
      {
        // Safe to expand
        if ( is_input( tfo ) )
          remove_input( tfo );
        make_tfo( tfo );
        window_.tfos.insert_unique( tfo );

        // Add fanouts as new TFOs
        ntk_.foreach_fanout( tfo, [&]( auto const& no ) {
          if ( !is_contained( no ) )
          {
            make_tfo( no );
            window_.tfos.insert_unique( no );
            front_nodes_.push_back( no );
          }
        } );

        // Promote new inputs
        for ( auto ni : new_inputs )
        {
          if ( !is_contained( ni ) )
          {
            make_input( ni );
            ntk_.foreach_output( ni, [&]( auto const& sig ) {
              window_.inputs.insert_unique( sig );
              window_.divs.insert_unique( sig );
            } );
          }
        }
      }

      // remove processed tfo from frontier
      for ( int i = (int)front_nodes_.size() - 1; i >= 0; --i )
      {
        if ( front_nodes_[i] == tfo )
        {
          front_nodes_.erase( i );
          break;
        }
      }
    }

    window_.outputs.clear();
    // --- 3) Determine outermost cut outputs ---
    for ( auto const& n : window_.tfos )
    {
      bool external_fanout = false;
      ntk_.foreach_fanout( n, [&]( auto const& fo ) {
        if ( !is_tfo( fo ) )
          external_fanout = true;
      } );
      if ( external_fanout || ntk_.is_po( n ) )
      {
        ntk_.foreach_output( n, [&]( auto const& sig ) {
          window_.outputs.insert_unique( sig );
        } );
      }
    }
    if ( window_.outputs.empty() )
    {
      ntk_.foreach_output( window_.pivot, [&]( auto const& f ) {
        window_.outputs.push_back( f );
      } );
    }
  }

  void tfi_expansion()
  {
    bool front_is_inputs = false;
    front_.clear();
    for ( auto const& l : window_.inputs )
      front_.push_back( l );

    while ( front_.size() <= ( ps_.max_num_leaves + ps_.max_reconvergence ) )
    {
      auto const leaf = choose_removable_leaf();
      if ( leaf )
      {
        make_contained( *leaf );
        ntk_.foreach_fanin( *leaf, [&]( auto const& fi, auto ii ) {
          auto const ni = ntk_.get_node( fi );
          if ( !is_contained( ni ) )
          {
            make_input( ni );
            ntk_.foreach_output( ni, [&]( auto const& f ) {
              front_.push_back( f );
            } );
          }
        } );

        for ( int i = front_.size() - 1; i >= 0; --i )
        {
          auto const n = ntk_.get_node( front_[i] );
          if ( !is_input( n ) )
            front_.erase( i );
        }

        if ( front_.size() <= ps_.max_num_leaves )
        {
          front_is_inputs = true;
          window_.inputs.clear();
          for ( auto const& l : front_ )
            window_.inputs.push_back( l );
        }
        else
        {
          front_is_inputs = false;
        }
      }
      else
        break;
    }

    for ( auto const& f : front_ )
    {
      auto const n = ntk_.get_node( f );
      make_alien( n );
    }
    for ( auto const& f : window_.inputs )
    {
      auto const n = ntk_.get_node( f );
      make_contained( n );
    }

    for ( auto const& f : window_.inputs )
    {
      auto const n = ntk_.get_node( f );
      if ( !is_input( n ) )
      {
        make_input( n );
        ntk_.foreach_output( n, [&]( auto const& fo ) {
          window_.divs.push_back( fo );
        } );
      }
    }

    if ( !front_is_inputs )
    {
      for ( auto const& l : window_.inputs )
      {
        auto const n = ntk_.get_node( l );
        ntk_.foreach_fanin( n, [&]( auto const& fi, auto i ) {
          auto const ni = ntk_.get_node( fi );
          cleanup_to_front( ni );
        } );
      }
    }

    ntk_.foreach_fanin( window_.pivot, [&]( auto const& fi, auto ii ) {
      auto const ni = ntk_.get_node( fi );
      if ( !is_input( ni ) )
        label_to_inputs( ni );
    } );
  }

  void side_expansion()
  {
    int start = 0;
    while ( start < window_.divs.size() )
    {
      int num_divs = window_.divs.size();
      for ( int i = start; i < num_divs; ++i )
      {
        auto const& n = ntk_.get_node( window_.divs[i] );
        ntk_.foreach_fanout( n, [&]( auto const& no ) {
          bool add = !is_contained( no );
          if ( add )
          {
            ntk_.foreach_fanin( no, [&]( auto const& fi, auto ii ) {
              auto const ni = ntk_.get_node( fi );
              add &= ( is_contained( ni ) && ( !is_mffc( ni ) ) && ( !is_tfo( ni ) ) );
            } );
          }
          if ( add )
          {
            make_divisor( no );
            ntk_.foreach_output( no, [&]( auto const& fo ) {
              window_.divs.push_back( fo );
            } );
          }
        } );
      }
      start = num_divs;
    }
  }

  bool is_input_dominated( node_index_t const& n )
  {
    bool dominated = true;
    ntk_.foreach_fanin( n, [&]( auto const& f, auto i ) {
      dominated &= is_contained( ntk_.get_node( f ) );
    } );
    return dominated;
  }

  void cleanup_to_front( node_index_t const& n )
  {
    if ( is_input( n ) )
      return;
    if ( !is_contained( n ) )
      return;
    if ( is_contained( n ) )
      make_alien( n );

    ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
      auto const ni = ntk_.get_node( fi );
      cleanup_to_front( ni );
    } );
  }

  void label_to_inputs( node_index_t const& n )
  {
    if ( !is_contained( n ) )
      return;
    if ( is_input( n ) )
      return;
    if ( is_mffc( n ) )
      return;
    if ( is_divisor( n ) )
      return;
    if ( in_mffc( n ) )
    {
      window_.mffc.push_back( n );
      make_mffc( n );
    }
    else
    {
      make_divisor( n );
      ntk_.foreach_output( n, [&]( auto const& fo ) {
        window_.divs.push_back( fo );
      } );
    }
    ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
      auto const ni = ntk_.get_node( fi );
      label_to_inputs( ni );
    } );
  }

  std::optional<node_index_t> choose_removable_leaf()
  {
    int best_cost = std::numeric_limits<int>::max();
    std::optional<node_index_t> best_node;
    for ( uint32_t i = 0; i < front_.size(); ++i )
    {
      auto n = ntk_.get_node( front_[i] );
      int cost = compute_leaf_cost( n );
      if ( cost < best_cost )
      {
        best_cost = cost;
        best_node = std::make_optional( n );
      }
    }
    return best_node;
  }

  std::optional<node_index_t> choose_removable_tfo()
  {
    int best_cost = std::numeric_limits<int>::max();
    bool best_promotes = false;
    std::optional<node_index_t> best_node;

    for ( uint32_t i = 0; i < front_nodes_.size(); ++i )
    {
      auto n = front_nodes_[i];
      auto c = compute_tfo_cost( n );

      if ( !best_node )
      {
        best_node = n;
        best_cost = c.delta_inputs;
        best_promotes = c.promotes_input;
        continue;
      }

      if ( c.promotes_input && !best_promotes )
      {
        // prioritize promotions over non-promotions
        best_node = n;
        best_cost = c.delta_inputs;
        best_promotes = true;
      }
      else if ( c.promotes_input == best_promotes && c.delta_inputs < best_cost )
      {
        // same promotion status, pick smaller delta
        best_node = n;
        best_cost = c.delta_inputs;
        best_promotes = c.promotes_input;
      }
    }

    return best_node;
  }

  bool in_mffc( node_index_t const& n )
  {
    if ( ntk_.is_pi( n ) || ntk_.is_po( n ) )
      return false;

    bool in = true;
    ntk_.foreach_fanout( n, [&]( auto const& no ) {
      in &= is_mffc( no );
    } );
    return in;
  }

  // cost to replace the node with its fanins
  int compute_leaf_cost( node_index_t const& n )
  {
    if ( ntk_.is_pi( n ) )
      return std::numeric_limits<int>::max();

    int cost = -static_cast<int>( ntk_.num_outputs( n ) );
    ntk_.foreach_fanin( n, [&]( auto const& fi, auto ii ) {
      auto const ni = ntk_.get_node( fi );
      if ( !is_contained( ni ) )
        cost += ntk_.num_outputs( ni );
    } );
    return cost;
  }

  struct tfo_cost_t
  {
    int delta_inputs;    // net input change
    bool promotes_input; // true if this expansion upgrades an input→TFO
  };

  tfo_cost_t compute_tfo_cost( node_index_t const& n )
  {
    std::set<node_index_t> new_inputs;

    ntk_.foreach_fanout( n, [&]( auto const& no ) {
      ntk_.foreach_fanin( no, [&]( auto const& fi, auto ) {
        auto const ni = ntk_.get_node( fi );
        if ( !is_contained( ni ) )
        {
          new_inputs.insert( ni );
        }
      } );
    } );

    int delta = static_cast<int>( new_inputs.size() );
    bool promote = false;

    if ( is_input( n ) )
    {
      delta -= 1;     // remove this input, add its fanins
      promote = true; // this is a reconvergence opportunity
    }

    return { delta, promote };
  }

private:
  void clear()
  {
    color_map_.resize();
    window_.clear();
  }

  void init( node_index_t const& n )
  {
    color_++;
    window_.pivot = n;
    // initialize the mffc with the pivot node
    window_.mffc.push_back( window_.pivot );
    make_mffc( window_.pivot );
    // initialize the outputs with the pivot signals
    ntk_.foreach_output( window_.pivot, [&]( auto const& f ) {
      window_.outputs.push_back( f );
    } );
    // initialize the leaves with the pivot's fanins
    ntk_.foreach_fanin( window_.pivot, [&]( auto const& fi, auto i ) {
      (void)i;
      auto const ni = ntk_.get_node( fi );
      if ( !is_input( ni ) )
      {
        make_input( ni );
        ntk_.foreach_output( ni, [&]( auto const& fo ) {
          window_.inputs.push_back( fo );
        } );
      }
    } );
  }

public:
  auto const& get_divisors() const
  {
    return window_.divs;
  }

  signal_t const& get_divisor( uint32_t const& index ) const
  {
    return window_.divs[index];
  }

  auto const& get_tfos() const
  {
    return window_.tfos;
  }

  auto const& get_outputs() const
  {
    return window_.outputs;
  }

  auto const& get_inputs() const
  {
    return window_.inputs;
  }

  auto const& get_divs() const
  {
    return window_.divs;
  }

  auto const& get_mffc() const
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

  template<uint32_t Cap>
  void topological_sort( rinox::math::static_vector<node_index_t, Cap>& nodes )
  {
    topo_++;
    static constexpr uint32_t maskTag = 0xFFFFF000u;   // upper 20 bits
    static constexpr uint32_t maskCount = 0x00000FFFu; // lower 12 bits
    static constexpr unsigned shiftTag = 12;
    uint32_t tag = topo_;

    // initialize indegree with tag and 0 count
    for ( auto const& n : nodes )
    {
      ntk_.set_value( n, ( tag << shiftTag ) | 0u );
    }

    // compute indegree counts (within this subset)
    for ( auto const& n : nodes )
    {
      ntk_.foreach_fanin( n, [&]( auto const& fi, auto ) {
        auto const& ni = ntk_.get_node( fi );
        if ( ( ( ntk_.value( ni ) >> shiftTag ) == tag ) )
        {
          uint32_t cnt = ( ntk_.value( n ) & maskCount ) + 1;
          assert( cnt <= maskCount && "indegree overflow (>4095 fanins)" );
          ntk_.set_value( n, ( tag << shiftTag ) | cnt );
        }
      } );
    }

    std::queue<node_index_t> q;
    for ( auto const& n : nodes )
    {
      if ( ( ntk_.value( n ) & maskCount ) == 0 )
        q.push( n );
    }

    rinox::math::static_vector<node_index_t, Cap> result;

    while ( !q.empty() )
    {
      auto n = q.front();
      q.pop();
      result.push_back( n );

      ntk_.foreach_fanout( n, [&]( auto const& no ) {
        if ( ( ntk_.value( no ) >> shiftTag ) == tag )
        {
          uint32_t cnt = ( ntk_.value( no ) & maskCount );
          assert( cnt > 0 && "indegree underflow in topological_sort" );
          cnt--;
          ntk_.set_value( no, ( tag << shiftTag ) | cnt );
          if ( cnt == 0 )
            q.push( no );
        }
      } );
    }

    if ( result.size() != nodes.size() )
    {
      throw std::runtime_error( "Cycle detected in topological_sort" );
    }

    nodes = std::move( result );
  }

  template<uint32_t Cap>
  void topological_sort( rinox::math::static_vector<signal_t, Cap>& signals )
  {
    rinox::math::static_vector<node_index_t, Cap> nodes;

    for ( auto const& s : signals )
    {
      if ( s.output == 0 )
        nodes.push_back( ntk_.get_node( s ) );
    }

    topological_sort<Cap>( nodes );

    rinox::math::static_vector<signal_t, Cap> result;
    for ( auto const& n : nodes )
    {
      ntk_.foreach_output( n, [&]( auto const& f ) {
        result.push_back( f );
      } );
    }

    signals = std::move( result );
  }

  void remove_redundant()
  {
    static constexpr uint32_t mffc_id = 1u;
    static constexpr uint32_t divs_id = 2u;
    static constexpr uint32_t inps_id = 3u;
    static constexpr uint32_t tfos_id = 4u;

    // mffc
    front_nodes_.clear();
    for ( int i = window_.mffc.size() - 1; i >= 0; --i )
    {
      auto const n = window_.mffc[i];
      if ( is_mffc( n ) )
      {
        front_nodes_.push_back( n );
        make_alien( n );
      }
    }
    window_.mffc.clear();
    for ( auto const& n : front_nodes_ )
    {
      window_.mffc.push_back( n );
      make_mffc( n );
    }
    // tfo
    front_nodes_.clear();
    for ( int i = window_.tfos.size() - 1; i >= 0; --i )
    {
      auto const n = window_.tfos[i];
      if ( is_tfo( n ) )
      {
        front_nodes_.push_back( n );
        make_alien( n );
      }
    }
    window_.tfos.clear();
    for ( auto const& n : front_nodes_ )
    {
      window_.tfos.push_back( n );
      make_tfo( n );
    }

    // inputs
    front_nodes_.clear();
    for ( int i = window_.inputs.size() - 1; i >= 0; --i )
    {
      auto const n = ntk_.get_node( window_.inputs[i] );
      if ( is_input( n ) )
      {
        front_nodes_.push_back( n );
        make_alien( n );
      }
    }
    window_.inputs.clear();
    for ( auto const& n : front_nodes_ )
    {
      ntk_.foreach_output( n, [&]( auto const& f ) {
        window_.inputs.push_back( f );
      } );
      make_input( n );
    }
    // divisors
    front_nodes_.clear();
    for ( int i = window_.divs.size() - 1; i >= 0; --i )
    {
      auto const n = ntk_.get_node( window_.divs[i] );
      if ( is_input( n ) )
      {
        front_nodes_.push_back( n );
        color_map_[n] = inps_id;
      }
      else if ( is_divisor( n ) )
      {
        front_nodes_.push_back( n );
        color_map_[n] = divs_id;
      }
    }
    window_.divs.clear();
    for ( auto const& n : front_nodes_ )
    {
      ntk_.foreach_output( n, [&]( auto const& f ) {
        window_.divs.push_back( f );
      } );
      color_map_[n] = ( color_map_[n] | ( color_ << 3u ) );
    }
  }

public:
  bool is_tfo( node_index_t const& n ) const
  {
    return is_contained( n ) && ( color_map_[n] == ( 4u | ( color_ << 3u ) ) );
  }

  bool is_contained( node_index_t const& n ) const
  {
    return ( n < color_map_.size() ) && color_map_.has( n ) && ( ( color_map_[n] >> 3u ) == color_ );
  }

  void make_contained( node_index_t const& n )
  {
    color_map_[n] = ( color_ << 3u );
  }
  void make_alien( node_index_t const& n )
  {
    color_map_[n] = 0;
  }

  uint32_t get_type( node_index_t const& n )
  {
    static constexpr uint32_t mask = 7u;
    return ( color_map_[n] & mask );
  }

  void make_mffc( node_index_t const& n )
  {
    color_map_[n] = ( 1u | ( color_ << 3u ) );
  }

  void make_divisor( node_index_t const& n )
  {
    color_map_[n] = ( 2u | ( color_ << 3u ) );
  }

  void make_input( node_index_t const& n )
  {
    color_map_[n] = ( 3u | ( color_ << 3u ) );
  }

  void make_tfo( node_index_t const& n )
  {
    color_map_[n] = ( 4u | ( color_ << 3u ) );
  }

  void make_output( node_index_t const& n )
  {
    color_map_[n] = ( 5u | ( color_ << 3u ) );
  }

  void make_unknown( node_index_t const& n )
  {
    color_map_[n] = ( 7u | ( color_ << 3u ) );
  }

  bool is_unknown( node_index_t const& n ) const
  {
    return is_contained( n ) && ( color_map_[n] == ( ( color_ << 3u ) | 7u ) );
  }

  void mark_contained()
  {
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

  bool is_input( node_index_t const& n ) const
  {
    return is_contained( n ) && ( color_map_[n] == ( 3u | ( color_ << 3u ) ) );
  }

  bool is_mffc( node_index_t const& n ) const
  {
    return is_contained( n ) && ( color_map_[n] == ( 1u | ( color_ << 3u ) ) );
  }

  bool is_divisor( node_index_t const& n ) const
  {
    return is_contained( n ) && ( color_map_[n] == ( 2u | ( color_ << 3u ) ) );
  }

  bool is_valid() const
  {
    return st_.valid;
  }
#pragma endregion

  size_t num_inputs() const
  {
    return window_.inputs.size();
  }

  size_t num_mffc() const
  {
    return window_.mffc.size();
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
        ntk_.foreach_fanin( d, [&]( auto const& fi, auto ii ) {
          if ( ii < ( ntk_.fanin_size( d.index ) - 1 ) )
            std::cout << "signal{ n[" << fi.index << "], " << fi.output << " } , ";
          else
            std::cout << "signal{ n[" << fi.index << "], " << fi.output << " } } , { ";
        } );
        for ( auto i = 0u; i < ( ids.size() - 1 ); ++i )
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
      ntk_.foreach_fanin( m, [&]( auto const& fi, auto ii ) {
        if ( ii < ( ntk_.fanin_size( m ) - 1 ) )
          std::cout << "signal{ n[" << fi.index << "], " << fi.output << " } , ";
        else
          std::cout << "signal{ n[" << fi.index << "], " << fi.output << " } } , { ";
      } );
      for ( auto i = 0u; i < ( ids.size() - 1 ); ++i )
      {
        std::cout << ids[i] << ", ";
      }
      std::cout << ids.back() << " } ).index;" << std::endl;
    }

    for ( auto const& m : window_.tfos )
    {
      auto const ids = ntk_.get_binding_ids( m );
      std::cout << "n[" << m << "] = ntk.create_node( ";
      ntk_.foreach_fanin( m, [&]( auto const& fi, auto ii ) {
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
  std::vector<signal_t> best_inputs_;
  front_signals_t front_;
  front_nodes_t front_nodes_;
  mockturtle::incomplete_node_map<uint32_t, Ntk> color_map_;
  uint32_t color_ = 1u;
  uint32_t topo_ = 1u;
  Params const& ps_;
  window_manager_stats& st_;
};

} // namespace windowing

} // namespace rinox
