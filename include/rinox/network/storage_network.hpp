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
  \file network_storage.hpp
  \brief Storage for network network specializing the operations on the nodes.
  \details
  This file defines the storage for the network network, which is a specialized
  data structure designed to handle multiple-output gates and their bindings.
  It includes methods for creating primary inputs and outputs, managing nodes,
  and handling the functional properties of the network. The encapsulation of
  the storage allows for efficient manipulation of the network while maintaining
  the flexibility to support various gate functionalities and bindings.
  \note This storage is designed to work with the `network_network` class, which
  provides a higher-level interface for interacting with the network.
  \note The `network_network` class uses this storage to manage the nodes, inputs,
  outputs, and the library of gates. It provides methods for creating nodes,
  replacing nodes, and querying the network's structure and functionality.

  \ingroup network_storage
  \see rinox::network::storage_node
  \see rinox::network::storage_signal
  \see rinox::network::storage_types

  \author Andrea Costamagna
*/

#pragma once

#include "../evaluation/chains.hpp"
#include "../libraries/augmented_library.hpp"
#include "storage_node.hpp"
#include "storage_signal.hpp"
#include "utils.hpp"
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/detail/foreach.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <queue>
#include <vector>

namespace rinox
{

/*!
 * \namespace network
 * \brief Types and utilities related to the network network data structure.
 */
namespace network
{

/*! \brief Compact storage for nodes in the network network.
 *
 * This structure represents the storage in network networks, enabling the
 * encapsulation of the detailed operations on nodes, inputs, and outputs.
 * It provides methods for creating primary inputs and outputs, managing nodes,
 * and handling the functional properties of the network. The storage is designed
 * to efficiently manage the nodes and their relationships, allowing for operations
 * such as creating nodes, replacing nodes, and querying the network's structure.
 *
 * \tparam NumBitsOutputs Number of bits used to represent the output pin specifier.
 */
template<design_type_t DesignType, uint32_t NumBitsOutputs>
class storage_network
{
#pragma region Types and constructors

public:
  using augmented_library_t = libraries::augmented_library<DesignType>;
  using gate = typename augmented_library_t::gate;
  using gate_t = typename augmented_library_t::gate_t;
  using list_t = evaluation::chains::large_xag_chain;
  using node_t = storage_node<NumBitsOutputs>;
  using signal_t = storage_signal<NumBitsOutputs>;

  /*! \brief The storage constructor.
   *
   * This constructor initializes the storage with a given library of gates.
   * It reserves space for a maximum number of nodes and initializes the first
   * two nodes as constants (0 and 1).
   *
   * \param gates The gates to be included in the library.
   */
  template<
      network::design_type_t D = DesignType,
      std::enable_if_t<D == network::design_type_t::CELL_BASED, int> = 0>
  storage_network( augmented_library_t library )
      : library( library )
  {
    /* reserve space for nodes */
    nodes.reserve( 10000u ); // reserve to avoid frequent reallocations

    /* we reserve the first two nodes for constants */
    init();
  }

  /*! \brief The storage constructor.
   *
   * This constructor initializes the storage with a given library of gates.
   * It reserves space for a maximum number of nodes and initializes the first
   * two nodes as constants (0 and 1).
   *
   * \param gates The gates to be included in the library.
   */
  template<
      network::design_type_t D = DesignType,
      std::enable_if_t<D == network::design_type_t::CELL_BASED, int> = 0>
  storage_network( std::vector<gate> gates )
      : library( gates )
  {
    /* reserve space for nodes */
    nodes.reserve( 10000u ); // reserve to avoid frequent reallocations

    /* we reserve the first two nodes for constants */
    init();
  }

  /*! \brief The storage constructor.
   *
   * This constructor initializes the storage with a given library of gates.
   * It reserves space for a maximum number of nodes and initializes the first
   * two nodes as constants (0 and 1).
   *
   * \param gates The gates to be included in the library.
   */
  storage_network()
  {
    /* reserve space for nodes */
    nodes.reserve( 10000u ); // reserve to avoid frequent reallocations

    /* we reserve the first two nodes for constants */
    init();
  }

  void init()
  {
    kitty::dynamic_truth_table const0( 0u );
    kitty::dynamic_truth_table const1( 0u );
    kitty::create_from_binary_string( const0, "0" );
    kitty::create_from_binary_string( const1, "1" );
    auto const id0 = library.get_id( const0 );
    auto const id1 = library.get_id( const1 );

    if ( id0 && id1 )
    {
      node_t new_node;
      new_node.children.clear();
      new_node.outputs = decltype( new_node.outputs )( 1u );
      auto const pin0 = output_pin_t{ *id0, pin_type_t::CONSTANT };
      auto const pin1 = output_pin_t{ *id1, pin_type_t::CONSTANT };
      new_node.outputs[0] = pin0;
      nodes.emplace_back( new_node );
      new_node.outputs[0] = pin1;
      nodes.emplace_back( new_node );
    }
    else
    {
      nodes.emplace_back();
      nodes.emplace_back();
    }
  }

#pragma endregion

#pragma region Primary I / O and constants

  /*! \brief Creates a constant signal.
   *
   * This method creates a signal representing a constant value (0 or 1).
   * It returns a signal with the appropriate index and output pin.
   *
   * \param value The constant value to be represented (true for 1, false for 0).
   * \return A signal representing the constant value.
   */
  signal_t get_constant( bool value ) const
  {
    return value ? signal_t{ 1, 0 } : signal_t{ 0, 0 };
  }

  /*! \brief Creates a primary input signal.
   *
   * This method creates a primary input signal and adds it to the storage.
   * It returns a signal with the index of the newly created primary input.
   * A PI stores its index in the only fanin it has.
   *
   * \return A signal representing the primary input.
   */
  signal_t create_pi()
  {
    const auto index = get_new_index();
    node_t input( pin_type_t::PI );
    input.children = { inputs.size() };
    nodes[index] = input;
    inputs.emplace_back( index );
    return signal_t{ index, 0 };
  }

  /*! \brief Creates a primary output signal.
   *
   * This method creates a primary output signal from a given signal.
   * It increases the reference count for the node to avoid incorrect deletions
   * and updates the output pin type to indicate it is a primary output. A node
   * can be used as PO more than once, so the number of fanouts is equal to the
   * fanout size of the output pins, plus the number of times one of its pins
   * is used as a primary output.
   * It returns the index of the newly created primary output.
   *
   * \param f The signal representing the function to be used as primary output.
   * \return The index of the newly created primary output.
   */
  uint32_t create_po( signal_t const& f )
  {
    /* increase ref-count to children */
    nodes[f.index].fanout_count++;
    nodes[f.index].outputs[f.output].type |= pin_type_t::PO;
    auto const po_index = static_cast<uint32_t>( outputs.size() );
    outputs.emplace_back( f.index, f.output );
    return po_index;
  }

  bool constant( node_index_t const& n ) const
  {
    return library.constant( get_binding_ids( n ) );
  }

  /*! \brief Check if the node is a constant */
  bool is_constant( node_index_t const& n ) const
  {
    auto const& pins = nodes[n].outputs;
    bool const inter = has_intersection( pins[0].type, pin_type_t::CONSTANT );
    bool const exter = inter || ( !is_pi( n ) && library.is_constant( get_binding_ids( n ) ) );
    return inter || exter;
  }

  /*! \brief Check if the node is a combinational input (CI)
   *
   * \param n The node to check.
   * \return True if the node is a combinational input, false otherwise.
   */
  bool is_ci( node_index_t const& n ) const
  {
    auto const& pins = nodes[n].outputs;
    return has_intersection( pins[0].type, pin_type_t::PI ) ||
           has_intersection( pins[0].type, pin_type_t::CI );
  }

  /*! \brief Check if the node is a primary input (PI)
   *
   * \param n The node to check.
   * \return True if the node is a primary input, false otherwise.
   */
  bool is_pi( node_index_t const& n ) const
  {
    return is_ci( n );
  }

  /*! \brief Check if the node is a primary output (PO)
   *
   * \param n The node to check.
   * \param output The output pin index to check (default is 0).
   * \return True if the node is a primary output, false otherwise.
   */
  bool is_po( node_index_t const& n, uint32_t output = 0 ) const
  {
    auto& pins = nodes[n].outputs;
    return has_intersection( pins[output].type, pin_type_t::PO );
  }

  bool is_po( signal_t f ) const
  {
    return is_po( f.index, f.output );
  }

  /*! \brief Returns if node index is not 0 */
  bool constant_value( node_index_t const& n ) const
  {
    return n != 0;
  }

#pragma endregion

#pragma region Create special functions
  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = true>
  signal_t create_not( signal_t const& a )
  {
    static uint64_t _not = 0x1;
    kitty::dynamic_truth_table tt_not( 1 );
    kitty::create_from_words( tt_not, &_not, &_not + 1 );
    return create_node<SignalsVec, DoStrash>( std::vector<signal_t>{ a }, tt_not );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = true>
  signal_t create_and( signal_t const& a, signal_t const& b )
  {
    static uint64_t _and = 0x8;
    kitty::dynamic_truth_table tt_and( 2 );
    kitty::create_from_words( tt_and, &_and, &_and + 1 );
    return create_node<SignalsVec, DoStrash>( std::vector<signal_t>{ a, b }, tt_and );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = true>
  signal_t create_nand( signal_t const& a, signal_t const& b )
  {
    static uint64_t _nand = 0xE;
    kitty::dynamic_truth_table tt_nand( 2 );
    kitty::create_from_words( tt_nand, &_nand, &_nand + 1 );
    return create_node<SignalsVec, DoStrash>( std::vector<signal_t>{ a, b }, tt_nand );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = true>
  signal_t create_or( signal_t const& a, signal_t const& b )
  {
    static uint64_t _or = 0xe;
    kitty::dynamic_truth_table tt_or( 2 );
    kitty::create_from_words( tt_or, &_or, &_or + 1 );
    return create_node<SignalsVec, DoStrash>( std::vector<signal_t>{ a, b }, tt_or );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = true>
  signal_t create_xor( signal_t const& a, signal_t const& b )
  {
    static uint64_t _xor = 0x6;
    kitty::dynamic_truth_table tt_xor( 2 );
    kitty::create_from_words( tt_xor, &_xor, &_xor + 1 );
    return create_node<SignalsVec, DoStrash>( std::vector<signal_t>{ a, b }, tt_xor );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = true>
  signal_t create_maj( signal_t const& a, signal_t const& b, signal_t const& c )
  {
    static uint64_t _maj = 0xe8;
    kitty::dynamic_truth_table tt_maj( 3 );
    kitty::create_from_words( tt_maj, &_maj, &_maj + 1 );
    return create_node<SignalsVec, DoStrash>( std::vector<signal_t>{ a, b, c }, tt_maj );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = true>
  signal_t create_ite( signal_t const& a, signal_t const& b, signal_t const& c )
  {
    static uint64_t _ite = 0xd8;
    kitty::dynamic_truth_table tt_ite( 3 );
    kitty::create_from_words( tt_ite, &_ite, &_ite + 1 );
    return create_node<SignalsVec, DoStrash>( std::vector<signal_t>{ a, b, c }, tt_ite );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = true>
  signal_t create_xor3( signal_t const& a, signal_t const& b, signal_t const& c )
  {
    static uint64_t _xor3 = 0x96;
    kitty::dynamic_truth_table tt_xor3( 3 );
    kitty::create_from_words( tt_xor3, &_xor3, &_xor3 + 1 );
    return create_node<SignalsVec, DoStrash>( std::vector<signal_t>{ a, b, c }, tt_xor3 );
  }

  unsigned int get_id( kitty::dynamic_truth_table const& tt ) const
  {
    auto id = library.get_id( tt );
    return *id;
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t create_node( SignalsVec const& children, kitty::dynamic_truth_table const& tt )
  {
    auto id = library.get_id( tt );
    if ( !id )
    {
      std::cerr << "[e] No binding found in the library.\n";
    }
    node_t n = create_storage_node( children, std::vector<unsigned int>{ *id } );
    return create_node( children, n );
  }
#pragma endregion

#pragma region Create arbitrary functions

  /*! \brief Create a detailed node to be stored.
   *
   * \param children input signals.
   * \param ids binding identifier of the output pins.
   */
  template<typename SignalsVec>
  node_t create_storage_node( SignalsVec const& children, std::vector<unsigned int> const& ids )
  {
    if constexpr ( DesignType == design_type_t::CELL_BASED )
    {
      for ( auto i = 1; i < ids.size(); ++i )
      {
        assert( library.get_name( ids[i] ) == library.get_name( ids[0] ) &&
                "Multiple-output nodes are expected to have the same name" );
      }
    }

    node_t new_node;
    std::copy( children.begin(), children.end(), std::back_inserter( new_node.children ) );

    new_node.outputs = decltype( new_node.outputs )( ids.size() );

    for ( auto i = 0; i < ids.size(); ++i )
      new_node.outputs[i] = { ids[i], pin_type_t::INTERNAL };

    return new_node;
  }

  /*! \brief Create a new node with multiple outputs.
   *
   * This method creates a new node with the specified children and output IDs.
   * It updates the fanout counts of the children and returns a signal representing
   * the new node.
   *
   * \param children The child signals that this node will depend on.
   * \param ids The IDs of the outputs for this node.
   * \return A signal representing the newly created node.
   */
  template<typename SignalsVec = std::vector<signal_t>>
  signal_t create_node( SignalsVec const& children, node_t const& n )
  {
    const auto index = get_new_index();
    nodes[index] = n;

    /* increase ref-count to children */
    for ( auto c : children )
    {
      nodes[c.index].fanout_count++;
      nodes[c.index].outputs[c.output].fanout_count++;
      nodes[c.index].outputs[c.output].fanout.push_back( index );
    }

    if ( hash.find( n ) != hash.end() )
    {
      hash[n].push_back( index );
    }
    else
    {
      hash[n] = { index };
    }

    return { index, 0 };
  }

  /*! \brief Create a new node from a truth table.
   *
   * This method creates a new node from its functionality, assuming LUT-like
   * values of the area and delays.
   *
   * \param children The child signals that this node will depend on.
   * \param function The truth table of the node's functionality.
   * \return The binding identifier of the gate stored in the library.
   */
  template<
      network::design_type_t D = DesignType,
      std::enable_if_t<D == network::design_type_t::ARRAY_BASED, int> = 0>
  uint32_t insert( kitty::dynamic_truth_table const& function )
  {
    return library.add_gate( function );
  }

  /*! \brief Create a new node from a vector of truth tables.
   *
   * This method creates a new node from its functionality, assuming LUT-like
   * values of the area and delays.
   *
   * \param children The child signals that this node will depend on.
   * \param functions The truth table of the node's functionality.
   * \return The binding identifier of the gate stored in the library.
   */
  template<
      network::design_type_t D = DesignType,
      std::enable_if_t<D == network::design_type_t::ARRAY_BASED, int> = 0>
  std::vector<uint32_t> insert( std::vector<kitty::dynamic_truth_table> const& functions )
  {
    return library.add_gates( functions );
  }

#pragma endregion

#pragma region Restructuring

  /*! \brief Update the list of POs when a signal in the old list is replaced.
   *
   * \param old_node The index of the old node to be replaced.
   * \param new_signals The new signals to replace the old node's outputs.
   */
  void replace_in_outputs( node_index_t const& old_node,
                           std::vector<signal_t> const& new_signals )
  {
    foreach_output_pin( old_node, [&]( auto const& pin, uint32_t i ) {
      signal_t const old_signal = signal_t{ old_node, i };
      if ( is_po( old_signal ) )
      {
        /* replace the output signals with the new signal */
        replace_output( old_signal, new_signals[i] );
      }
    } );
  }

  /*! \brief Replace an output signal in the outputs list.
   *
   * This method replaces an old signal with a new signal in the outputs list.
   * It updates the fanout count of the new signal and the old signal,
   *
   * \param old_signal The old signal to be replaced.
   * \param new_signal The new signal to replace the old signal.
   */
  void replace_output( signal_t const& old_signal,
                       signal_t const& new_signal )
  {
    bool found = false;
    for ( auto& output : outputs )
    {
      if ( output == old_signal )
      {
        found = true;
        /* replace the old signal with the new signal in the outputs */
        output = new_signal;
        /* update the fanout count of the new signal */
        nodes[new_signal.index].fanout_count++;
        nodes[old_signal.index].fanout_count--;
        nodes[old_signal.index].outputs[old_signal.output].type &= ~pin_type_t::PO;
        nodes[new_signal.index].outputs[new_signal.output].type |= pin_type_t::PO;
      }
    }
    assert( found && "Output signal not found in the outputs list" );
  }

  /*! \brief Insert a fanout node for a signal.
   *
   * This method inserts a fanout for a given signal and node index.
   * It updates the fanout count of the signal and the output pin accordingly.
   *
   * \param f The signal for which the fanout is to be inserted.
   * \param n The node index to be added as a fanout.
   */
  void insert_fanout( signal_t const& f, node_index_t const& n )
  {
    auto& fanout = nodes[f.index].outputs[f.output].fanout;
    uint32_t const occurrences = std::count( fanout.begin(),
                                             fanout.end(),
                                             n );
    if ( occurrences > 0 )
    {
      /* if the fanout already exists, we do not need to insert it again */
      return;
    }
    nodes[f.index].fanout_count++;
    nodes[f.index].outputs[f.output].fanout_count++;
    nodes[f.index].outputs[f.output].fanout.push_back( n );
  }

  /*! \brief Delete a fanout node for a signal.
   *
   * This method deletes a fanout for a given signal and node index.
   * It updates the fanout count of the signal and the output pin accordingly.
   *
   * \param f The signal for which the fanout is to be deleted.
   * \param n The node index to be removed from the fanout.
   */
  void delete_fanout( signal_t const& f, node_index_t const& n )
  {
    auto& fanout = nodes[f.index].outputs[f.output].fanout;
    uint32_t const occurrences = std::count( fanout.begin(),
                                             fanout.end(),
                                             n );
    nodes[f.index].fanout_count -= occurrences;
    nodes[f.index].outputs[f.output].fanout_count -= occurrences;
    fanout.erase( std::remove( fanout.begin(),
                               fanout.end(),
                               n ),
                  fanout.end() );
  }

  /*! \brief Update the interconnections of a node.
   *
   * This method updates the fanin-fanout information in the network by replacing an old signal
   * with a new signal in the fanin of a specified node. It updates the fanout count of the new signal
   * and applies events to notify about the modification.
   *
   * \param root The index of the node where the replacement occurs.
   * \param old_signal The old signal to be replaced.
   * \param new_signal The new signal to replace the old node.
   */
  void update_nets( node_index_t const& root,
                    signal_t const& old_signal,
                    signal_t new_signal )
  {
    auto& nd_root = nodes[root];
    for ( auto& child : nd_root.children )
    {
      if ( child == old_signal )
      {
        /* add the root node to the new signal's fanout */
        insert_fanout( new_signal, root );

        /* update the fanout of the old signal */
        delete_fanout( child, root );

        /* replace the old signal with the new signal in the children */
        child = new_signal;
      }
    }
  }

  /*! \brief Delete a node from the network.
   *
   * This method removes a node from the network, marking it as dead and
   * updating the fanout counts of its children. It also triggers events
   * for the deletion of the node.
   *
   * \param n The index of the node to be removed.
   */
  void delete_node( node_index_t const& n )
  {
    /* remove the node from the hash table if present */
    auto it = hash.find( nodes[n] );
    if ( it != hash.end() )
    {
      auto& list = it->second;
      list.erase( std::remove( list.begin(), list.end(), n ), list.end() );
      if ( list.empty() ) // commented
      {
        hash.erase( it );
      }
    }
    /* mark the node as dead */
    for ( auto& pin : nodes[n].outputs )
    {
      pin.type |= pin_type_t::DEAD;
      pin.fanout.clear();
      pin.fanout_count = 0;
    }
    nodes[n].fanout_count = 0;
    nodes[n].children.clear();

    dead_nodes.push( n );
  }

#pragma endregion

#pragma region Structural properties

  /*! \brief Returns true since the network is combinational.
   *
   * TODO: Add support for sequential elements in the future.
   * \return Always returns true.
   */
  bool is_combinational() const
  {
    return true;
  }

  bool is_multioutput( node_index_t const& n ) const
  {
    return num_outputs( n ) > 1;
  }

  template<bool DoStrash = false,
           network::design_type_t D = DesignType,
           std::enable_if_t<D == network::design_type_t::CELL_BASED, int> = 0>
  bool is_multioutput( std::string const& name ) const
  {
    return library.is_multioutput( name );
  }

  /*! \brief Check if the node is dead
   *
   * \param n The node to check.
   * \return True if the node is dead, false otherwise.
   *
   * A dead node is one where all output pins are marked as DEAD.
   * This typically indicates that the node is no longer used in the network.
   */
  bool is_dead( node_index_t const& n ) const
  {
    bool all_dead{ true };
    bool one_dead{ false };
    for ( auto const& pin : nodes[n].outputs )
    {
      bool const dead = has_intersection( pin.type, pin_type_t::DEAD );
      all_dead &= dead;
      one_dead |= dead;
    }
    assert( !( all_dead ^ one_dead ) );
    /* A dead node is simply a dangling node */
    return all_dead;
  }

  inline bool is_constant( signal_t const& f ) const
  {
    auto const& outputs = nodes[f.index].outputs;
    uint32_t i = f.output;
    return outputs.size() > 0 && has_intersection( outputs[i].type, pin_type_t::CONSTANT );
  }

  auto size() const
  {
    return static_cast<uint32_t>( nodes.size() );
  }

  auto num_cis() const
  {
    return static_cast<uint32_t>( inputs.size() );
  }

  auto num_cos() const
  {
    return static_cast<uint32_t>( outputs.size() );
  }

  auto num_pis() const
  {
    return static_cast<uint32_t>( inputs.size() );
  }

  auto num_pos() const
  {
    return static_cast<uint32_t>( outputs.size() );
  }

  auto num_gates() const
  {
    return static_cast<uint32_t>( nodes.size() - inputs.size() - 2 );
  }

  uint32_t num_outputs( node_index_t const& n ) const
  {
    return static_cast<uint32_t>( nodes[n].outputs.size() );
  }

  uint32_t fanin_size( node_index_t const& n ) const
  {
    return static_cast<uint32_t>( nodes[n].children.size() );
  }

  uint32_t fanout_size( node_index_t const& n ) const
  {
    return nodes[n].fanout_count;
  }

  uint32_t incr_fanout_size( node_index_t const& n )
  {
    return nodes[n].fanout_count++;
  }

  uint32_t decr_fanout_size( node_index_t const& n )
  {
    return --nodes[n].fanout_count;
  }

  uint32_t incr_fanout_size_pin( node_index_t const& n, uint32_t pin_index )
  {
    return ++nodes[n].outputs[pin_index].fanout_count;
  }

  uint32_t decr_fanout_size_pin( node_index_t const& n, uint32_t pin_index )
  {
    return --nodes[n].outputs[pin_index].fanout_count;
  }

  uint32_t fanout_size_pin( node_index_t const& n, uint32_t pin_index ) const
  {
    return nodes[n].outputs[pin_index].fanout_count;
  }

  bool is_function( node_index_t const& n ) const
  {
    auto const& outputs = nodes[n].outputs;
    return ( outputs.size() > 0 ) && ( has_intersection( outputs[0].type, network::pin_type_t::INTERNAL ) ||
                                       has_intersection( outputs[0].type, network::pin_type_t::PO ) );
  }

  /*! \brief Checks if a given node is already present in the storage.
   *
   * Uses structural hashing to check if a node is in the network.
   */
  std::optional<node_index_t> find( node_t const& n ) const
  {
    auto it = hash.find( n );
    if ( it != hash.end() )
    {
      assert( !is_dead( it->second[0] ) &&
              "The node should not be dead when looking for it" );
      return it->second[0];
    }
    return std::nullopt;
  }

  /*! \brief Checks if a node is in the fanin of another one.
   */
  bool in_fanin( node_index_t parent, node_index_t other ) const
  {
    bool in_fanin = false;
    auto& nobj = nodes[parent];
    for ( auto& f : nobj.children )
    {
      if ( f.index == other )
      {
        in_fanin = true;
        return in_fanin;
      }
    }

    return in_fanin;
  }

#pragma endregion

#pragma region Functional properties
  kitty::dynamic_truth_table signal_function( const signal_t& f ) const
  {
    auto const& outputs = nodes[f.index].outputs;
    auto const& id = outputs[f.output].id;
    return library[id].function;
  }
#pragma endregion

#pragma region Nodes and signals

  node_index_t ci_at( uint32_t index ) const
  {
    assert( index < inputs.size() );
    return *( inputs.begin() + index );
  }

  signal_t co_at( uint32_t index ) const
  {
    assert( index < outputs.size() );
    return *( outputs.begin() + index );
  }

  node_index_t pi_at( uint32_t index ) const
  {
    assert( index < inputs.size() );
    return *( inputs.begin() + index );
  }

  signal_t po_at( uint32_t index ) const
  {
    assert( index < outputs.size() );
    return *( outputs.begin() + index );
  }

  uint32_t pi_index( node_index_t const& n ) const
  {
    auto const& outputs = nodes[n].outputs;
    assert( has_intersection( outputs[0].type, network::pin_type_t::PI ) );
    return static_cast<uint32_t>( nodes[n].children[0].data );
  }

  uint32_t po_index( signal_t const& f ) const
  {
    for ( uint32_t i = 0; i < outputs.size(); ++i )
    {
      if ( f == outputs[i] )
        return i;
    }
    assert( false && "PO does not exist" );
    return std::numeric_limits<uint32_t>::max();
  }
#pragma endregion

#pragma region Node and signal iterators

  template<typename Fn>
  void foreach_node( Fn&& fn ) const
  {
    for ( node_index_t n = 2u; n < nodes.size(); ++n )
    {
      if ( !is_dead( n ) )
      {
        fn( n );
      }
    }
  }

  template<typename Fn>
  void foreach_ci( Fn&& fn ) const
  {
    mockturtle::detail::foreach_element( inputs.begin(), inputs.end(), fn );
  }

  template<typename Fn>
  void foreach_co( Fn&& fn ) const
  {
    mockturtle::detail::foreach_element( outputs.begin(), outputs.end(), fn );
  }

  template<typename Fn>
  void foreach_pi( Fn&& fn ) const
  {
    mockturtle::detail::foreach_element( inputs.begin(), inputs.end(), fn );
  }

  template<typename Fn>
  void foreach_po( Fn&& fn ) const
  {
    mockturtle::detail::foreach_element( outputs.begin(), outputs.end(), fn );
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    auto const nodes_size = nodes.size();
    auto r = mockturtle::range<uint64_t>( 2u, nodes_size ); /* start from 2 to avoid constants */
    mockturtle::detail::foreach_element_if(
        r.begin(), r.end(),
        [this]( auto n ) { return !is_ci( n ) && !is_dead( n ); },
        fn );
  }

  template<typename Fn>
  void foreach_fanin( node_index_t const& n, Fn&& fn ) const
  {
    if ( n <= 1 || is_ci( n ) )
      return;

    mockturtle::detail::foreach_element( nodes[n].children.begin(), nodes[n].children.end(), fn );
  }

  /*! \brief Iterate over the output pins of a node.
   *
   * This method iterates over the output pins of a specified node and applies
   * a function to each pin. The function receives the pin and its index as arguments.
   *
   * \param n The index of the node whose output pins are to be iterated.
   * \param fn The function to apply to each output pin.
   */
  template<typename Fn>
  void foreach_output_pin( node_index_t const& n, Fn&& fn ) const
  {
    auto& pins = nodes[n].outputs;
    for ( auto i = 0u; i < pins.size(); ++i )
    {
      fn( pins[i], i );
    }
  }

  /*! \brief Iterate over the outputs of a node.
   *
   * This method iterates over the output pins of a specified node represented
   * as signals, applying a function to each signal.
   *
   * \param n The index of the node whose output pins are to be iterated.
   * \param fn The function to apply to each output pin.
   */
  template<typename Fn>
  void foreach_output( node_index_t const& n, Fn&& fn ) const
  {
    auto& pins = nodes[n].outputs;
    for ( auto i = 0u; i < pins.size(); ++i )
    {
      auto const f = signal_t{ n, i };
      fn( f );
    }
  }

  template<typename Fn>
  void foreach_fanout( output_pin_t const& pin, Fn&& fn ) const
  {
    auto const fanout = pin.fanout;
    for ( uint32_t i = 0; i < fanout.size(); ++i )
    {
      node_index_t const& n = fanout[i];
      fn( n, i );
    }
  }

  template<typename Fn>
  void foreach_fanout( signal_t const& f, Fn&& fn ) const
  {
    auto& fanout = nodes[f.index].outputs[f.output].fanout;
    for ( uint32_t i = 0; i < fanout.size(); ++i )
    {
      node_index_t const& n = fanout[i];
      fn( n );
    }
  }

  template<typename Fn>
  void foreach_fanout( node_index_t const& n, Fn&& fn ) const
  {
    foreach_output_pin( n, [&]( auto const& pin, auto i ) {
      foreach_fanout( pin, [&]( auto const& fanout_node, auto j ) {
        fn( fanout_node );
      } );
    } );
  }

#pragma endregion

#pragma region Custom node values
  void clear_values()
  {
    std::for_each( nodes.begin(), nodes.end(), []( auto& n ) { n.user_data = 0; } );
  }

  uint32_t value( node_index_t const& n ) const
  {
    return nodes[n].user_data;
  }

  void set_value( node_index_t const& n, uint32_t v )
  {
    nodes[n].user_data = v;
  }

  void set_name( signal_t const& f, std::string const& name )
  {
    names_map[f] = name;
  }

  uint32_t incr_value( node_index_t const& n )
  {
    return static_cast<uint32_t>( nodes[n].user_data++ );
  }

  uint32_t decr_value( node_index_t const& n )
  {
    return static_cast<uint32_t>( --nodes[n].user_data );
  }

  bool has_output_name( uint32_t po_index )
  {
    return output_names.size() > po_index;
  }

  std::string get_output_name( uint32_t po_index )
  {
    return output_names[po_index];
  }

  bool has_input_name( uint32_t pi_index )
  {
    return names_map.find( make_signal( pi_at( pi_index ) ) ) != names_map.end();
  }

  std::string get_input_name( uint32_t pi_index )
  {
    signal_t const f{ pi_at( pi_index ), 0 };
    return get_name( f );
  }

  void set_output_name( uint32_t po_index, std::string const& name )
  {
    output_names.resize( outputs.size() );
    output_names[po_index] = name;
  }

  void set_input_name( uint32_t pi_index, std::string const& name )
  {
    signal_t const f{ pi_at( pi_index ), 0 };
    names_map[f.data] = name;
  }

  bool has_name( signal_t const& f )
  {
    if ( names_map.find( f.data ) != names_map.end() )
      return true;
    if ( is_po( f ) )
      return output_names.size() == num_pos();
    return false;
  }

  std::string get_name( signal_t const& f )
  {
    return names_map[f];
  }

  std::string get_network_name() const
  {
    return module_name;
  }

  void set_network_name( std::string const& name )
  {
    module_name = name;
  }

#pragma endregion

#pragma region Visited flags
  void clear_visited()
  {
    std::for_each( nodes.begin(), nodes.end(), []( auto& n ) { n.traversal_id = 0; } );
  }

  auto visited( node_index_t const& n ) const
  {
    return nodes[n].traversal_id;
  }

  void set_visited( node_index_t const& n, uint32_t v )
  {
    nodes[n].traversal_id = v;
  }

  uint32_t get_trav_id() const
  {
    return trav_id;
  }

  unsigned int get_pin_id( std::string const& gate_name, std::string const& pin_name ) const
  {
    return library.get_pin_id( gate_name, pin_name );
  }

  void incr_trav_id()
  {
    if ( trav_id > ( std::numeric_limits<uint32_t>::max() - 10 ) )
    {
      std::cout << "[w] Traversal identifier exceeded safe treshold. Forced reset" << std::endl;
      clear_values();
      clear_visited();
      trav_id = 0;
    }
    ++trav_id;
  }
#pragma endregion

#pragma region Getters
  /*! \brief Get the children of a node.
   *
   * This method retrieves the children of a specified node.
   * It returns a vector of signals representing the children.
   *
   * \param n The index of the node whose children are to be retrieved.
   * \return A vector of signals representing the children of the node.
   */
  std::vector<signal_t> const& get_children( node_index_t const& n ) const
  {
    return nodes[n].children;
  }

  /*! \brief Get the fanin of a node.
   *
   * This method retrieves the fanin of a specified node.
   * It returns a vector of node indeces in the immediate fanin.
   *
   * \param n The index of the node whose children are to be retrieved.
   * \return A vector of node indices representing the children of the node.
   */
  std::vector<node_index_t> const get_fanins( node_index_t const& n ) const
  {
    std::vector<node_index_t> fanin;
    for ( auto const& child : nodes[n].children )
    {
      fanin.push_back( child.index );
    }
    return fanin;
  }

  /*! \brief
   */
  list_t const& get_chain( uint32_t id ) const
  {
    return library.get_chain( id );
  }

  double const& get_area( node_index_t const& n ) const
  {
    auto const& g = get_binding( signal_t{ n, 0 } );
    return g.area;
  }

  std::vector<kitty::dynamic_truth_table> get_functions( node_index_t const& n ) const
  {
    std::vector<kitty::dynamic_truth_table> tts;
    foreach_output( n, [&]( auto const& f ) {
      auto const& g = get_binding( signal_t{ n, 0 } );
      tts.push_back( g.function );
    } );
    return tts;
  }

  std::string const get_expression( signal_t const& f ) const
  {
    auto const& g = get_binding( f );
    return g.expression;
  }

  /*! \brief Get the binding identifiers of the output pins in a node
   */
  std::vector<uint32_t> get_binding_ids( node_index_t const& n ) const
  {
    std::vector<uint32_t> ids;
    for ( auto const& pin : nodes[n].outputs )
      ids.push_back( pin.id );
    return ids;
  }

  /*! \brief Get the binding identifiers of the output pins in a node
   */
  std::vector<unsigned int> get_binding_ids( std::string const& name ) const
  {
    return library.get_binding_ids( name );
  }

  auto const& get_binding( signal_t const& f ) const
  {
    auto const& pin = nodes[f.index].outputs[f.output];
    return library.get_gate( pin.id );
  }

  double get_max_pin_delay( signal_t const& f, uint32_t i ) const
  {
    auto const& pin = nodes[f.index].outputs[f.output];
    return library.get_max_pin_delay( pin.id, i );
  }

  double get_min_pin_delay( signal_t const& f, uint32_t i ) const
  {
    auto const& pin = nodes[f.index].outputs[f.output];
    return library.get_min_pin_delay( pin.id, i );
  }

  double get_input_load( signal_t const& f, uint32_t i ) const
  {
    auto const& pin = nodes[f.index].outputs[f.output];
    return library.get_input_load( pin.id, i );
  }

  std::vector<gate_t> const& get_library() const
  {
    return library.get_aug_gates();
  }

  node_index_t get_new_index()
  {
    // if ( dead_nodes.empty() )
    {
      nodes.emplace_back( pin_type_t::NONE );
      return static_cast<node_index_t>( nodes.size() - 1 );
    }
    // auto n = dead_nodes.front();
    // dead_nodes.pop();
    // nodes[n].clear();
    // nodes[n] = node_t( pin_type_t::NONE );
    // return n;
  }

  uint32_t get_fanin_number( unsigned int id, std::string const& pin_name ) const
  {
    return library.get_fanin_number( id, pin_name );
  }
#pragma endregion

#pragma region Bindings

  bool has_binding( signal_t const& f ) const
  {
    return has_binding( f.index );
  }

  bool has_binding( node_index_t const& n ) const
  {
    return ( !is_ci( n ) && !is_pi( n ) );
  }

  bool has_gate( std::string const& name ) const
  {
    return library.has_gate( name );
  }

  bool is_input_pin( std::string const& gate_name, std::string const& pin_name ) const
  {
    return library.is_input_pin( gate_name, pin_name );
  }

  bool is_output_pin( std::string const& gate_name, std::string const& pin_name ) const
  {
    return library.is_output_pin( gate_name, pin_name );
  }
#pragma endregion

  /*! \brief Traversal ID for graph algorithms.
   *
   * This ID is used to mark nodes during traversal operations.
   * It is initialized to zero and can be incremented for each traversal.
   */
  uint32_t trav_id = 0u;

  /*! \brief The nodes in the network network.
   *
   * This vector stores all the nodes in the network network, each represented
   * by a `storage_node` object. It includes primary inputs, outputs, and
   * internal nodes.
   */
  std::vector<node_t> nodes;

  /*! \brief The nodes that were killed in the network network */
  std::queue<node_index_t> dead_nodes;

  /*! \brief The primary inputs of the network network.
   *
   * This vector stores the indices of the primary input nodes in the network.
   * Each input corresponds to a node that can be used as a starting point for
   * logic operations.
   */
  std::vector<node_index_t> inputs;

  /*! \brief The primary outputs of the network network.
   *
   * This vector stores the signals representing the primary outputs of the network.
   * Each output corresponds to a signal that can be used to observe the results
   * of logic operations in the network.
   */
  std::vector<signal_t> outputs;
  std::vector<std::string> output_names;
  std::unordered_map<uint64_t, std::string> names_map;

  /*! \brief The library of gates used in the network network.
   *
   * This library contains the gates that can be used to create nodes in the network.
   * It is initialized with a set of gates and provides methods for accessing and
   * manipulating the gates, as well as an AIG list representation for simulation.
   */
  augmented_library_t library;
  std::string module_name;

  /*! \brief Hash map for fast node lookups.
   *
   * This hash map allows for quick access to nodes based on their indices.
   * It uses a custom hash function to ensure efficient storage and retrieval
   * of nodes in the network.
   */
  phmap::flat_hash_map<node_t, std::vector<node_index_t>, network::node_hash<NumBitsOutputs>> hash;
};

} // namespace network

} // namespace rinox
