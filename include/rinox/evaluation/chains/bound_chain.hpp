/* rinox: C++ logic network library
 * Copyright (C) 2025  EPFL
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
  \file bound_chain.hpp
  \brief List of indices to represent small Majority Inverter Graphs (MIG)

  \author Andrea Costamagna
*/

#pragma once

#include "../../boolean/boolean.hpp"
#include "../../math/math.hpp"
#include "../../network/utils.hpp"
#include <cassert>
#include <mockturtle/networks/block.hpp>
#include <vector>

namespace rinox
{

namespace evaluation
{

namespace chains
{

/*! \brief Boolean chain of gates from a technology library.
 *
 * The inputs are associated with the literals 0, ..., num_inputs - 1.
 * The subsequent literals identify the nodes in the chain.
 */
template<network::design_type_t DesignType>
class bound_chain
{
public:
  using element_type = uint32_t;

private:
  /*! \brief Node of a mapped index chain */
  struct node_t
  {
    node_t() = default;
    node_t( const node_t& ) = default;
    node_t( node_t&& ) noexcept = default;
    node_t& operator=( const node_t& ) = default;
    node_t& operator=( node_t&& ) noexcept = default;

    node_t( std::vector<element_type> const& fanins, uint32_t id )
        : fanins( fanins ), id( id )
    {}

    bool operator==( node_t const& other ) const
    {
      return fanins == other.fanins && id == other.id;
    }

    bool operator!=( node_t const& other ) const
    {
      return fanins != other.fanins || id == other.id;
    }

    std::vector<element_type> fanins; /**< literals of fanins */
    uint32_t id;                      /**< binding id */
  };

public:
#pragma region Types and Constructors

  bound_chain( uint32_t num_inputs, size_t reserve_size )
      : num_inputs( num_inputs )
  {
    nodes.reserve( reserve_size );
  }

  explicit bound_chain( uint32_t num_inputs )
      : bound_chain( num_inputs, 10u )
  {}

  bound_chain()
      : bound_chain( 0 )
  {}

  ~bound_chain() = default;
  bound_chain( bound_chain const& ) = default;
  bound_chain( bound_chain&& ) noexcept = default;
  bound_chain& operator=( bound_chain const& ) = default;
  bound_chain& operator=( bound_chain&& ) noexcept = default;

  void clear()
  {
    nodes.clear();
    outputs.clear();
  }

#pragma endregion

#pragma region Equality

  bool operator==( bound_chain const& other ) const
  {
    if ( num_inputs != other.num_inputs || outputs != other.outputs || nodes.size() != other.nodes.size() )
      return false;

    for ( size_t i = 0; i < nodes.size(); ++i )
    {
      if ( nodes[i] != other.nodes[i] )
        return false;
    }
    return true;
  }

#pragma endregion

#pragma region Primary I/O and node creation

  void add_inputs( uint32_t n = 1 )
  {
    num_inputs += n;
  }

  void add_output( element_type const v )
  {
    outputs.push_back( v );
  }

  element_type pi_at( uint32_t index ) const
  {
    assert( index < num_inputs );
    return index;
  }

  element_type po_at( uint32_t index ) const
  {
    assert( index < outputs.size() );
    return outputs[index];
  }

  bool is_pi( element_type f ) const
  {
    return f < num_inputs;
  }

  /*! \brief Create a node in the chain. The returned literal uniquely identifies it. */
  element_type add_gate( std::vector<element_type> const& fanins, uint32_t id )
  {
    const element_type f = static_cast<element_type>( nodes.size() + num_inputs );
    nodes.emplace_back( fanins, id );
    return f;
  }

  /*! \brief Replace in node */
  void replace_in_node( uint32_t node, uint32_t fanin, element_type other )
  {
    auto& n = nodes[node];
    n.fanins[fanin] = other;
  }

  /*! \brief Replace in node */
  void replace_output( uint32_t index, element_type other )
  {
    outputs[index] = other;
  }
#pragma endregion

#pragma region Iterators

  template<typename Fn>
  void foreach_pi( Fn&& fn ) const
  {
    for ( uint32_t i = 0; i < num_inputs; ++i )
    {
      fn( static_cast<element_type>( i ) );
    }
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    for ( size_t i = 0; i < nodes.size(); ++i )
    {
      auto const& n = nodes[i];
      fn( n.fanins, n.id, i );
    }
  }

  template<typename Fn>
  void foreach_gate_rev( Fn&& fn ) const
  {
    for ( int i = static_cast<int>( nodes.size() ) - 1; i >= 0; --i )
    {
      auto const& n = nodes[i];
      fn( n.fanins, n.id, i );
    }
  }

  template<typename Fn>
  void foreach_po( Fn&& fn ) const
  {
    for ( size_t i = 0; i < outputs.size(); ++i )
    {
      fn( outputs[i], i );
    }
  }

#pragma endregion

#pragma region Structural properties

  [[nodiscard]] uint32_t num_gates() const
  {
    return static_cast<uint32_t>( nodes.size() );
  }

  [[nodiscard]] uint32_t num_pis() const
  {
    return num_inputs;
  }

  [[nodiscard]] uint32_t num_pos() const
  {
    return static_cast<uint32_t>( outputs.size() );
  }

  [[nodiscard]] uint32_t size() const
  {
    return num_inputs + static_cast<uint32_t>( nodes.size() );
  }

#pragma endregion

#pragma region Getters

  template<typename Lib>
  [[nodiscard]] double get_area( Lib const& lib ) const
  {
    double area = 0;
    foreach_gate( [&]( auto const& fanin, auto const& id, auto i ) {
      area += lib.get_area( id );
    } );
    return area;
  }

  [[nodiscard]] std::vector<node_t> get_nodes() const
  {
    return nodes;
  }

  [[nodiscard]] std::vector<element_type> get_outputs() const
  {
    return outputs;
  }

  [[nodiscard]] element_type get_num_inputs() const
  {
    return num_inputs;
  }

  [[nodiscard]] auto get_pi_index( element_type lit ) const
  {
    return static_cast<uint32_t>( lit );
  }

  [[nodiscard]] auto get_node_index( element_type lit ) const
  {
    return static_cast<uint32_t>( lit ) - num_inputs;
  }

#pragma endregion

private:
  std::vector<node_t> nodes;
  std::vector<element_type> outputs;
  element_type num_inputs = 0;
};

/*! \brief Returns the longest path from the inputs to any output */
template<network::design_type_t DesignType, typename Lib>
std::vector<double> get_longest_paths( bound_chain<DesignType>& chain, Lib& library )
{
  std::vector<double> input_delays( chain.num_pis(), 0 );
  std::vector<double> nodes_delays( chain.num_pis() + chain.num_gates(), 0 );
  chain.foreach_gate_rev( [&]( auto const& fanin, auto id, auto i ) {
    for ( auto j = 0u; j < fanin.size(); ++j )
    {
      nodes_delays[fanin[j]] = std::min( nodes_delays[fanin[j]], nodes_delays[chain.num_pis() + i] - library.get_max_pin_delay( id, j ) );
      if ( chain.is_pi( fanin[j] ) )
        input_delays[fanin[j]] = nodes_delays[fanin[j]];
    }
  } );
  for ( auto& d : input_delays )
    d = -d;
  return input_delays;
}

/*! \brief Permutes the input variables to have the ones closest to the output first
 *
 *    *
 *   ***
 *   \**
 *    \**
 *   0 \*
 *    1 \*
 *     2
 *      3
 *
 * */
template<network::design_type_t DesignType, typename Lib>
void time_canonize( bound_chain<DesignType>& chain, Lib& library, boolean::symmetries_t const& symm )
{
  std::vector<uint8_t> inputs( chain.num_pis() );
  std::iota( inputs.begin(), inputs.end(), 0 );

  auto delays = get_longest_paths( chain, library );

  for ( uint8_t i = 0; i < chain.num_pis(); ++i )
  {
    if ( symm.has_symmetries( i ) )
    {
      uint8_t k = i;
      int j = i - 1;
      double delay = delays[inputs[i]];
      bool swapped = true;
      while ( swapped && ( j >= 0 ) )
      {
        if ( symm.symmetric( inputs[k], inputs[j] ) )
        {
          if ( delay > delays[j] )
          {
            std::swap( inputs[k], inputs[j] );
            std::swap( delays[k], delays[j] );
            k = j;
            swapped = true;
          }
          else
          {
            swapped = false;
          }
        }
        j--;
      }
    }
  }

  boolean::permutation_t perm( inputs );

  chain.foreach_gate( [&]( auto& fanin, auto id, auto i ) {
    for ( auto j = 0u; j < fanin.size(); ++j )
    {
      auto const lit = fanin[j];
      if ( chain.is_pi( lit ) )
      {
        chain.replace_in_node( i, j, perm.inverse( lit ) );
      }
    }
  } );

  chain.foreach_po( [&]( auto& lit, auto i ) {
    if ( chain.is_pi( lit ) )
    {
      chain.replace_output( i, perm.inverse( lit ) );
    }
  } );
}

template<network::design_type_t DesignType>
void perm_canonize( bound_chain<DesignType>& chain, boolean::permutation_t const& perm )
{
  chain.foreach_gate( [&]( auto& fanins, auto const& id, auto i ) {
    for ( int j = 0; j < fanins.size(); ++j )
    {
      if ( chain.is_pi( fanins[j] ) )
      {
        auto k = perm.inverse( fanins[j] );
        chain.replace_in_node( i, j, k );
      }
    }
  } );

  chain.foreach_po( [&]( auto& output, auto i ) {
    if ( chain.is_pi( output ) )
    {
      auto k = perm.inverse( output );
      chain.replace_output( i, k );
    }
  } );
}

// leaves must be a signal container
template<typename Ntk, network::design_type_t DesignType, typename Signals, bool DoStrash = true>
typename Ntk::signal_t insert( Ntk& ntk,
                               Signals const& inputs,
                               bound_chain<DesignType> const& chain )
{
  using signals_t = std::vector<typename Ntk::signal_t>;
  signals_t fs;

  for ( auto i = 0; i < chain.num_pis(); ++i )
    fs.emplace_back( inputs[i] );

  chain.foreach_gate( [&]( auto const& fanin, auto id, auto i ) {
    signals_t children( fanin.size() );
    for ( auto i = 0u; i < fanin.size(); ++i )
    {
      children[i] = fs[fanin[i]];
    }
    fs.push_back( ntk.template create_node<signals_t, DoStrash>( children, id ) );
  } );
  return fs[chain.po_at( 0 )];
}

template<typename Ntk, network::design_type_t DesignType, typename Signals = std::vector<typename Ntk::signal>>
void extract( bound_chain<DesignType>& chain,
              Ntk& ntk,
              Signals const& inputs,
              typename Ntk::signal const& output )
{
  using signal = typename Ntk::signal;
  using node = typename Ntk::node;
  using value_type = typename bound_chain<DesignType>::element_type;

  if ( chain.num_pis() < inputs.size() )
  {
    std::cerr << "[e] not enough PIs in the chain" << std::endl;
    return;
  }

  std::unordered_map<uint64_t, value_type> sig_to_lit;

  // Recursive construction
  std::function<void( signal )> construct_rec = [&]( signal f ) {
    node n = ntk.get_node( f );

    if ( ntk.visited( n ) == ntk.trav_id() )
      return;

    if ( ntk.is_pi( n ) )
    {
      std::cerr << "[e] unexpected unmarked PI in logic cone\n";
      return;
    }

    std::vector<value_type> children;

    ntk.foreach_fanin( n, [&]( auto const& fi, auto ) {
      construct_rec( fi );
      children.emplace_back( sig_to_lit[fi] );
    } );

    value_type v;
    if constexpr ( mockturtle::has_has_cell_v<Ntk> )
    {
      auto cell = ntk.get_cell( n );
      assert( cell.gates.size() <= 1 && "[e] multiple output support not available" );
      auto const id = cell.gates[0].id;
      v = chain.add_gate( children, id );
      sig_to_lit[f] = v;
    }
    else
    {
      auto const ids = ntk.get_binding_ids( n );
      for ( auto i = 0u; i < ids.size(); ++i )
      {
        v = chain.add_gate( children, ids[i] );
        sig_to_lit[signal( { n, i } )] = v;
      }
    }

    ntk.set_visited( n, ntk.trav_id() );
  };

  // Assign values to input nodes
  ntk.incr_trav_id();
  for ( std::size_t i = 0; i < inputs.size(); ++i )
  {
    if ( inputs[i].data > std::numeric_limits<uint32_t>::max() )
      continue;
    auto const n = ntk.get_node( inputs[i] );
    ntk.set_visited( n, ntk.trav_id() );
    sig_to_lit[inputs[i]] = static_cast<value_type>( i ); // assign PI index
  }

  construct_rec( output );
  chain.add_output( sig_to_lit[output] );
}

} // namespace chains

} // namespace evaluation

} // namespace rinox
