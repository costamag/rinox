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
  \file bound.hpp
  \brief Bound network with multiple-output gates support.

  This data structure is a general logic representation which can be used to
  represent a network type with multiple-output cells. Natively, this network
  representation is designed for enabling efficient optimization after technology
  mapping. However, its generality allows us to use it to represent any network
  type in rinox, including AIGs, XAIGs, MIGs, XMGs, etc. To support these
  representations, signals can be complemented.

  \author Andrea Costamagna
*/

#pragma once

#include "../evaluation/evaluation.hpp"
#include "../libraries/libraries.hpp"
#include "signal_map.hpp"
#include "storage_network.hpp"
#include "storage_node.hpp"
#include "storage_signal.hpp"
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/detail/foreach.hpp>
#include <mockturtle/networks/events.hpp>
#include <mockturtle/networks/storage.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/algorithm.hpp>
#include <mockturtle/utils/truth_table_cache.hpp>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>

#include <algorithm>
#include <memory>

namespace rinox
{

namespace network
{

/*! \brief Network of gates from a technology library.
 *
 * \tparam MaxNumOutputs Maximum number of outputs of the cells in the library.
 */
template<design_type_t DesignType = design_type_t::CELL_BASED, uint32_t MaxNumOutputs = 2u>
class bound_network
{
public:
#pragma region Types and constructors

  static constexpr auto NumBitsOutputs = bits_required<MaxNumOutputs>();
  static constexpr design_type_t design_t = DesignType;
  static constexpr bool is_bound_network_type = true;

  /* aliases used in this class */
  using storage_t = std::shared_ptr<storage_network<DesignType, NumBitsOutputs>>;
  using list_t = evaluation::chains::large_xag_chain;
  using node_t = storage_node<NumBitsOutputs>;
  using signal_t = storage_signal<NumBitsOutputs>;
  using hash_t = signal_hash<NumBitsOutputs>;
  using node_index_t = uint64_t;

  /* aliases for compatibility with the other network types */
  static constexpr auto min_fanin_size = 1;
  static constexpr auto max_fanin_size = 32;
  static constexpr auto max_num_outputs = 1u << NumBitsOutputs;
  static constexpr auto kMaxNumOutputs = MaxNumOutputs;
  using base_type = bound_network<DesignType, MaxNumOutputs>;
  using storage = storage_t;
  using signal = signal_t;
  using node = node_index_t;

  /*! \brief Constructor from a technology library.
   *
   * \param gates The gates in the technology library.
   */
  template<
      design_type_t D = DesignType,
      std::enable_if_t<D == design_type_t::CELL_BASED, int> = 0>
  bound_network( libraries::augmented_library<DesignType> library )
      : _storage( std::make_shared<storage_network<DesignType, NumBitsOutputs>>( library ) ),
        _events( std::make_shared<typename decltype( _events )::element_type>() )
  {
  }

  /*! \brief Constructor from a technology library.
   *
   * \param gates The gates in the technology library.
   */
  template<
      design_type_t D = DesignType,
      std::enable_if_t<D == design_type_t::CELL_BASED, int> = 0>
  bound_network( std::vector<mockturtle::gate> gates )
      : _storage( std::make_shared<storage_network<DesignType, NumBitsOutputs>>( gates ) ),
        _events( std::make_shared<typename decltype( _events )::element_type>() )
  {
  }

  /*! \brief Constructor for array-based designs.
   *
   * \param gates The gates in the technology library.
   */
  template<
      design_type_t D = DesignType,
      std::enable_if_t<D == design_type_t::ARRAY_BASED, int> = 0>
  bound_network()
      : _storage( std::make_shared<storage_network<DesignType, NumBitsOutputs>>() ),
        _events( std::make_shared<typename decltype( _events )::element_type>() )
  {
  }

  /*! \brief Constructor from a storage object.
   *
   * This constructor is used to create a bound network from an existing storage
   * object, allowing for cloning and manipulation of the network without
   * needing to recreate the storage structure.
   *
   * \param storage The storage object containing the network data.
   */
  bound_network( std::shared_ptr<storage_network<DesignType, NumBitsOutputs>> storage )
      : _storage( storage ),
        _events( std::make_shared<typename decltype( _events )::element_type>() )
  {
  }

  /*! \brief Clone the current network.
   *
   * This method creates a new instance of the bound network with a copy of the
   * current storage. It is useful for creating a separate instance of the
   * network that can be modified independently of the original.
   *
   * \return A new bound_network instance with cloned storage.
   */
  bound_network<DesignType, MaxNumOutputs> clone() const
  {
    return { std::make_shared<storage_network<DesignType, NumBitsOutputs>>( *_storage ) };
  }

#pragma endregion

#pragma region Primary I / O and constants
public:
  /*! \brief Returns a constant signal.
   *
   * This method returns a signal representing a constant value (0 or 1).
   * The value can be specified as an argument, with the default being false (0).
   *
   * \param value The constant value to be represented (true for 1, false for 0).
   * \return A signal representing the constant value.
   */
  signal_t get_constant( bool value = false ) const
  {
    return _storage->get_constant( value );
  }

  /*! \brief Creates a primary input signal.
   *
   * \return A signal representing the primary input.
   */
  signal_t create_pi()
  {
    return _storage->create_pi();
  }

  /*! \brief Label a signal as primary output.
   *
   * \param f The signal to be added to the primary outputs.
   * \return A unique identifier for the primary output.
   */
  uint32_t create_po( signal_t const& f )
  {
    return _storage->create_po( f );
  }

  /*! \brief Check if a node is a constant.
   *
   * \param n The node index to check.
   * \return True if the node is a constant, false otherwise.
   */
  bool is_constant( node_index_t n ) const
  {
    return _storage->is_constant( n );
  }

  bool constant( node_index_t n ) const
  {
    return _storage->constant( n );
  }

  /*! \brief Check if a node is a combinational input (CI).
   *
   * \param n The node index to check.
   * \return True if the node is a combinational input, false otherwise.
   */
  bool is_ci( node_index_t n ) const
  {
    return _storage->is_ci( n );
  }

  /*! \brief Check if a node is a primary input (PI).
   *
   * \param n The node index to check.
   * \return True if the node is a primary input, false otherwise.
   */
  bool is_pi( node_index_t n ) const
  {
    assert(_storage && "storage unavailable");
    bool const res = _storage->is_pi( n );;
    return res;
  }

  bool is_inv(node_index_t n) const
  {
    if (is_pi(n) || (fanin_size(n)!=1)||(num_outputs(n)>1)) return false;
    auto const tt = node_function(n,0);
    return !kitty::is_normal(tt);
  }

  bool is_buf(node_index_t n) const
  {
    if (is_pi(n) || (fanin_size(n)!=1)||(num_outputs(n)>1)) return false;
    auto const tt = node_function(n,0);
    return kitty::is_normal(tt);
  }

  /*! \brief Check if a node is a primary output (PO).
   *
   * \param n The node index to check.
   * \param output The output pin index to check (default is 0).
   * \return True if the node is a primary output, false otherwise.
   */
  bool is_po( node_index_t n, uint32_t output ) const
  {
    return _storage->is_po( n, output );
  }

  bool is_po( node_index_t n ) const
  {
    bool po = false;
    foreach_output( n, [&]( auto const& f ) {
      po |= is_po( f );
    } );
    return po;
  }

  /*! \brief Check if a signal is a primary output (PO).
   *
   * \param f The signal to check.
   * \return True if the signal is a primary output, false otherwise.
   */
  bool is_po( signal_t const& f ) const
  {
    return is_po( f.index, f.output );
  }

  /*! \brief Check if a node is a constant 0 or not.
   *
   * \param n The node index to check.
   * \return False if the node represents the constant 0, true otherwise.
   */
  bool constant_value( node_index_t n ) const
  {
    return _storage->constant_value( n );
  }
#pragma endregion

#pragma region Create special functions
  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t create_not( signal_t const& a )
  {
    signal_t const fnew = _storage->template create_not<SignalsVec, DoStrash>( a );
    for ( auto const& fn : _events->on_add )
    {
      ( *fn )( fnew.index );
    }
    return fnew;
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t create_and( signal_t const& a, signal_t const& b )
  {
    return _storage->template create_and<SignalsVec, DoStrash>( a, b );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t create_nand( signal_t const& a, signal_t const& b )
  {
    return _storage->template create_nand<SignalsVec, DoStrash>( a, b );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t create_or( signal_t const& a, signal_t const& b )
  {
    return _storage->template create_or<SignalsVec, DoStrash>( a, b );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t create_xor( signal_t const& a, signal_t const& b )
  {
    return _storage->template create_xor<SignalsVec, DoStrash>( a, b );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t create_maj( signal_t const& a, signal_t const& b, signal_t const& c )
  {
    return _storage->template create_maj<SignalsVec, DoStrash>( a, b, c );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t create_ite( signal_t const& a, signal_t const& b, signal_t const& c )
  {
    return _storage->template create_ite<SignalsVec, DoStrash>( a, b, c );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t create_xor3( signal_t const& a, signal_t const& b, signal_t const& c )
  {
    return _storage->template create_xor3<SignalsVec, DoStrash>( a, b, c );
  }

#pragma endregion

#pragma region Create arbitrary functions

  /*! \brief Create a node from the fanin signals and binding IDs.
   *
   * This method creates a new node in the network with the specified children
   * and binding IDs. when more than one binding ID is provided, the node is a
   * multiple-output node, allowing for multiple outputs pins to be associated
   * with different functions.
   *
   * \param children The children signals of the new node.
   * \param ids The binding IDs for the outputs of the new node.
   * \tparam DoStrash If true, the node will be created with strashing enabled.
   * \return A signal representing the newly created node.
   */
  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t create_node( SignalsVec const& children,
                        std::vector<uint32_t> const& ids )
  {
    node_t const& n = _storage->create_storage_node( children, ids );

    /* structural hashing */
    if constexpr ( DoStrash )
    {
      const auto it = _storage->find( n );
      if ( it )
      {
        if ( !is_dead( *it ) )
        {
          return { *it, 0 };
        }
      }
    }

    signal_t const f = _storage->create_node( children, n );

    /* initialize the application specific value to 0 */
    set_value( f.index, 0 );

    for ( auto const& fn : _events->on_add )
    {
      ( *fn )( f.index );
    }
    return f;
  }

  /*! \brief Create a node from the truth table of its functionality.
   *
   * This method creates a new node from its functionality. If no function
   * exists in the database containing the desired functionality, it is inserted
   * in the truth table cache as a technology-independent gate for LUT-like
   * storage. This functionality allows bound networks to represent LUT networks
   *
   * \param children The children signals of the new node.
   * \param tt Truth table of the new node.
   * \tparam DoStrash If true, the node will be created with strashing enabled.
   * \return A signal representing the newly created node.
   */
  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false,
           design_type_t D = DesignType,
           std::enable_if_t<D == design_type_t::ARRAY_BASED, int> = 0>
  signal_t create_node( SignalsVec const& children,
                        kitty::dynamic_truth_table& tt )
  {
    auto const binding_id = _storage->insert( tt );
    return create_node<SignalsVec, DoStrash>( children, binding_id );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false,
           design_type_t D = DesignType,
           std::enable_if_t<D == design_type_t::ARRAY_BASED, int> = 0>
  signal_t create_node( SignalsVec const& children,
                        std::vector<kitty::dynamic_truth_table>& tts )
  {
    std::vector<unsigned int> ids;
    for ( auto const& tt : tts )
      ids.push_back(_storage->insert( tt ) );
    return create_node<SignalsVec, DoStrash>( children, ids );
  }

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false,
           design_type_t D = DesignType,
           std::enable_if_t<D == design_type_t::CELL_BASED, int> = 0>
  signal_t create_node( SignalsVec const& children,
                        kitty::dynamic_truth_table const& tt )
  {
    auto const binding_id = _storage->get_id( tt );
    if ( !binding_id )
      throw std::runtime_error( "[ERROR] function not found" );
    return create_node<SignalsVec, DoStrash>( children, binding_id );
  }

  unsigned int get_id( kitty::dynamic_truth_table const& tt ) const
  {
    return _storage->get_id( tt );
  }

  uint32_t insert( kitty::dynamic_truth_table const& function )
  {
    return _storage->insert( function );
  }

  /*! \brief Create a node from the truth table of its functionality.
   *
   * This method creates a new node from its functionality. If no function
   * exists in the database containing the desired functionality, it is inserted
   * in the truth table cache as a technology-independent gate for LUT-like
   * storage. This functionality allows bound networks to represent LUT networks
   *
   * \param children The children signals of the new node.
   * \param tt Truth table of the new node.
   * \tparam DoStrash If true, the node will be created with strashing enabled.
   * \return A signal representing the newly created node.
   */
  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false,
           design_type_t D = DesignType,
           std::enable_if_t<D == design_type_t::ARRAY_BASED, int> = 0>
  signal_t create_node( SignalsVec const& children,
                        std::vector<kitty::dynamic_truth_table> const& tts )
  {
    std::vector<uint32_t> const binding_ids = _storage->insert( tts );
    return create_node<SignalsVec, DoStrash>( children, binding_ids );
  }

  /*! \brief Create a node with a single binding ID.
   *
   * This method is a convenience overload for creating a node with a single
   * binding ID. It calls the more general `create_node` method with a vector
   * containing the single ID.
   *
   * \param children The children signals of the new node.
   * \param id The binding ID for the output of the new node.
   * \tparam DoStrash If true, the node will be created with strashing enabled.
   * \return A signal representing the newly created node.
   */
  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t create_node( std::vector<signal_t> const& children,
                        uint32_t id )
  {
    return create_node<SignalsVec, DoStrash>( children, std::vector<uint32_t>{ id } );
  }

  /*! \brief Clone a node from another bound network.
   *
   * This method creates a new node in the current network by cloning an existing
   * node from another bound network. It takes the source node and its children
   * signals, and creates a new node with the same binding IDs.
   *
   * \param other The other bound network from which to clone the node.
   * \param source The index of the source node to clone.
   * \param children The children signals of the new node.
   * \return A signal representing the newly cloned node.
   */

  template<typename SignalsVec = std::vector<signal_t>, bool DoStrash = false>
  signal_t clone_node( bound_network const& other,
                       node_index_t source,
                       std::vector<signal_t> const& children )
  {
    assert( !children.empty() );
    std::vector<uint32_t> const ids = other.get_binding_ids( source );
    return create_node<SignalsVec, DoStrash>( children, ids );
  }

#pragma endregion

#pragma region Restructuring

  /*! \brief Substitute a node with signals equivalent to its output pins.
   *
   * This method replaces a node's output pins with functionally equivalent
   * signals. It updates the fanout of the old node's outputs to point to the
   * new signals, effectively substituting the old node in the network.
   *
   * [ pin 0 ] -> [ new_signals[0] ]
   * ...
   * [ pin j ] -> [ new_signals[j] ]
   *
   * \param old_node The index of the old node to be replaced.
   * \param new_signals The new signals to replace the old node's outputs.
   */
  void substitute_node( node_index_t old_node,
                        std::vector<signal_t> const& new_signals )
  {
    /* update the signals to be used as primary outputs.
     * Highest priority so that on_modified events operate on the correct POs.
     */
    _storage->replace_in_outputs( old_node, new_signals );

    /* update the fanins/fanout information and trigger modified events */
    replace_in_node( old_node, new_signals );

    /* remove the node and trigger on_delete events */
    take_out_node( old_node, new_signals );
  }

  /*! \brief Substitute a node with a new signal.
   *
   * This method replaces an old node with a new signal in the network, falling
   * back to the general case of vectorized substitute node.
   *
   * \param old_node The index of the old node to be replaced.
   * \param new_signal The new signal to replace the old node's outputs.
   */
  void substitute_node( node_index_t old_node, signal_t const& new_signal )
  {
    /* fall back to the general case */
    substitute_node( old_node, std::vector<signal_t>{ new_signal } );
  }

  /*! \brief Update the fanin-fanout information in the network.
   *
   * Iterate over all output pins of the node to be removed and replace
   * the old node's outputs with the new signals. This method updates the
   * fanout count of the new signals and adjusts its fanout accordingly.
   *
   * \param old_node The index of the old node to be replaced.
   * \param new_signals The new signals to replace the old node's outputs.
   */
  void replace_in_node( node_index_t old_node,
                        std::vector<signal_t> const& new_signals )
  {
    auto old_outputs = num_outputs( old_node );
    if ( new_signals.size() != old_outputs )
    {
      throw std::runtime_error( "[ERROR] substitute_node: output count mismatch" );
    }

    _storage->foreach_output_pin( old_node, [&]( auto const& pin, auto i ) {
      signal_t const old_signal{ old_node, i };
      auto const& new_signal = new_signals.at( i );

      _storage->foreach_fanout( pin, [&]( auto const& fanout_node, auto j ) {
        (void)j;
        if ( new_signal.output >= num_outputs( new_signal.index ) )
        {
          throw std::out_of_range( "[ERROR] substitute_node: new_signal pin index out of range" );
        }
        replace_in_node( fanout_node, old_signal, new_signal );
      } );
    } );
  }

  /*! \brief Replace a node in the fanin of another node.
   *
   * This method replaces an old node with a new signal in the fanin of a specified node.
   * It updates the fanout count of the new signal and adjusts the outputs accordingly.
   *
   * \param root The index of the node where the replacement occurs.
   * \param old_signal The old signal to be replaced.
   * \param new_signal The new signal to replace the old node.
   */
  void replace_in_node( node_index_t root,
                        signal_t const& old_signal,
                        signal_t new_signal )
  {
    node_index_t const old_node = old_signal.index;

    if ( !_storage->in_fanin( root, old_node ) )
      return;

    auto const old_children = _storage->get_children( root );

    _storage->update_nets( root, old_signal, new_signal );

    /* provide the root and the old signals to the modified event.
     * This corresponds to all the information, since the new children can be
     * computed from the root.
     */
    for ( auto const& fn : _events->on_modified )
    {
      ( *fn )( root, old_children );
    }
  }

  /*! \brief Take out a node if it is not reused in the new nodes.
   *
   * This method checks if the old node is still used in the network.
   * If it is not, it removes the node and updates the fanout counts of its children.
   *
   * \param old_node The index of the old node to be removed.
   * \param new_signals The new signals to replace the old node's outputs.
   */
  void take_out_node( node_index_t old_node,
                      std::vector<signal_t> const& new_signals )
  {
    /* take out the node if it is not reused in the new nodes */
    for ( auto const& f : new_signals )
    {
      if ( f.index == old_node )
      {
        /* if the old node is still used, we cannot take it out */
        return;
      }
    }
    take_out_node( old_node );
  }

  void remove_edge(node_index_t from, node_index_t to)
  {
    _storage->remove_edge(from, to);
  }
  /*! \brief Take out a node from the network.
   *
   * This method removes a node from the network, marking it as dead and
   * updating the fanout counts of its children. It also triggers events
   * for the deletion of the node.
   *
   * \param n The index of the node to be removed.
   */
  void take_out_node( node_index_t n )
  {
    /* we cannot delete CIs, constants, or already dead nodes */
    if ( is_constant( n ) || is_ci( n ) || is_dead( n ) )
      return;

    auto children = _storage->get_children( n );

    /* NOTE: the node's information is not cleared-up yet, so we can
     * access the node's outputs or the node's fanins. Not its old fanouts.
     */
    for ( auto const& fn : _events->on_delete )
    {
      ( *fn )( n );
    }

    /* mark the node as dead */
    _storage->delete_node( n );

    /* if the node has been deleted, then deref fanout_size of
       fanins and try to take them out if their fanout_size become 0 */
    for ( auto i = 0; i < children.size(); ++i )
    {
      auto& child = children[i];
      _storage->delete_fanout( child, n );
      if ( fanout_size( child.index ) == 0 )
      {
        take_out_node( child.index );
      }
    }
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

  bool is_multioutput( node_index_t n ) const
  {
    return _storage->is_multioutput( n );
  }

  template<bool DoStrash = false,
           design_type_t D = DesignType,
           std::enable_if_t<D == design_type_t::CELL_BASED, int> = 0>
  bool is_multioutput( std::string const& name ) const
  {
    return _storage->is_multioutput( name );
  }

  inline bool is_dead( node_index_t n ) const
  {
    return _storage->is_dead( n );
  }

  auto size() const
  {
    return _storage->size();
  }

  auto signal_size() const
  {
    return size() << (NumBitsOutputs+2);
  }

  auto num_cis() const
  {
    return _storage->num_cis();
  }

  auto num_cos() const
  {
    return _storage->num_cos();
  }

  auto num_pis() const
  {
    return _storage->num_pis();
  }

  auto num_pos() const
  {
    return _storage->num_pos();
  }

  auto num_gates() const
  {
    return _storage->num_gates();
  }

  uint32_t num_outputs( node_index_t n ) const
  {
    return _storage->num_outputs( n );
  }

  uint32_t fanin_size( node_index_t n ) const
  {
    return _storage->fanin_size( n );
  }

  uint32_t fanout_size( node_index_t n ) const
  {
    return _storage->fanout_size( n );
  }

  uint32_t fanout_size( signal_t const& f ) const
  {
    return fanout_size_pin( f.index, f.output );
  }

  uint32_t incr_fanout_size( node_index_t n ) const
  {
    return _storage->incr_fanout_size( n );
  }

  uint32_t decr_fanout_size( node_index_t n ) const
  {
    return _storage->decr_fanout_size( n );
  }

  uint32_t incr_fanout_size_pin( node_index_t n, uint32_t pin_index ) const
  {
    return _storage->incr_fanout_size_pin( n, pin_index );
  }

  uint32_t decr_fanout_size_pin( node_index_t n, uint32_t pin_index ) const
  {
    return _storage->decr_fanout_size_pin( n, pin_index );
  }

  uint32_t fanout_size_pin( node_index_t n, uint32_t pin_index ) const
  {
    return _storage->fanout_size_pin( n, pin_index );
  }

  bool is_function( node_index_t n ) const
  {
    return _storage->is_function( n );
  }

  double area()
  {
    double total_area = 0;
    foreach_gate( [&]( auto const& n ) {
      if ( !is_dead( n ) )
        total_area += get_area( n );
    } );
    return total_area;
  }

  uint32_t CountGates() const
  {
    uint32_t count = 0;
    foreach_gate( [&]( auto const& n ){
      if ( !is_dead( n ) )
        count+=1;
    } );
    return count;
  }


  uint32_t CountDuals() const
  {
    uint32_t count = 0;
    foreach_gate( [&]( auto const& n ){
      if ( !is_dead( n ) && (num_outputs(n)==2) )
        count+=1;
    } );
    return count;
  }
#pragma endregion

#pragma region Functional properties
  kitty::dynamic_truth_table signal_function( const signal_t& f ) const
  {
    return _storage->signal_function( f );
  }

  kitty::dynamic_truth_table node_function( const node_index_t& n, uint32_t pin_index = 0 ) const
  {
    signal f = make_signal( n, pin_index );
    return signal_function( f );
  }
#pragma endregion

#pragma region Nodes and signals
  node_index_t get_node( signal_t const& f ) const
  {
    return f.index;
  }

  signal_t make_signal( node_index_t n, uint32_t output_pin ) const
  {
    return signal_t{ n, output_pin };
  }

  signal_t make_signal( node_index_t n ) const
  {
    return make_signal( n, 0 );
  }

  bool is_complemented( signal_t const& f ) const
  {
    (void)f;
    return false;
  }

  uint32_t get_output_pin( signal_t const& f ) const
  {
    return static_cast<uint32_t>( f.output );
  }

  signal_t next_output_pin( signal_t const& f ) const
  {
    return signal_t{ f.index, f.output + 1 };
  }

  uint32_t node_to_index( node const& n ) const
  {
    return static_cast<uint32_t>( n );
  }

  node_index_t index_to_node( uint32_t index ) const
  {
    return index;
  }

  uint64_t signal_to_index( signal_t const& f ) const
  {
    return static_cast<uint64_t>( f.data );
  }

  node_index_t ci_at( uint32_t index ) const
  {
    return _storage->ci_at( index );
  }

  signal_t co_at( uint32_t index ) const
  {
    return _storage->co_at( index );
  }

  node_index_t pi_at( uint32_t index ) const
  {
    return _storage->pi_at( index );
  }

  signal_t po_at( uint32_t index ) const
  {
    return _storage->po_at( index );
  }

  uint32_t pi_index( node_index_t n ) const
  {
    return _storage->pi_index( n );
  }

  uint32_t po_index( signal_t const& f ) const
  {
    return _storage->po_index( f );
  }
#pragma endregion

#pragma region Node and signal iterators
  template<typename Fn>
  void foreach_node( Fn&& fn ) const
  {
    _storage->foreach_node( fn );
  }

  template<typename Fn>
  void foreach_ci( Fn&& fn ) const
  {
    _storage->foreach_ci( fn );
  }

  template<typename Fn>
  void foreach_co( Fn&& fn ) const
  {
    _storage->foreach_co( fn );
  }

  template<typename Fn>
  void foreach_pi( Fn&& fn ) const
  {
    _storage->foreach_pi( fn );
  }

  template<typename Fn>
  void foreach_po( Fn&& fn ) const
  {
    _storage->foreach_po( fn );
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    _storage->foreach_gate( fn );
  }

  template<typename Fn>
  void foreach_fanin( node_index_t n, Fn&& fn ) const
  {
    _storage->foreach_fanin( n, fn );
  }

  template<typename Fn>
  void foreach_fanin( signal_t const& f, Fn&& fn ) const
  {
    auto const n = get_node( f );
    _storage->foreach_fanin( n, fn );
  }

  template<typename Fn>
  void foreach_fanout( node_index_t n, Fn&& fn ) const
  {
    _storage->foreach_fanout( n, fn );
  }

  template<typename Fn>
  void foreach_fanout( signal_t const& f, Fn&& fn ) const
  {
    _storage->foreach_fanout( f, fn );
  }

  template<typename Fn>
  void foreach_tfo_node( node_index_t n, Fn&& fn ) const
  {
    _storage->foreach_tfo_node( n, fn );
  }

  template<typename Fn>
  void foreach_output_pin( node_index_t n, Fn&& fn ) const
  {
    _storage->foreach_output_pin( n, fn );
  }

  template<typename Fn>
  void foreach_output( node_index_t n, Fn&& fn ) const
  {
    _storage->foreach_output( n, fn );
  }
#pragma endregion

#pragma region Simulate values

  /*! \brief Get the cached simulator for AIG index lists.
   *
   * Caching an unique simulator avoids reallocations of different simulation
   * engines, ensuring memory efficiency.
   */
  template<typename TT>
  std::shared_ptr<evaluation::chain_simulator<list_t, TT>> get_simulator() const
  {
    using simulator_t = evaluation::chain_simulator<list_t, TT>;
    static const std::shared_ptr<simulator_t> sim = std::make_shared<simulator_t>();
    return sim;
  }

  /*! \brief Simulation of the input patterns using the node's function.
   *
   * \param n index of the node to simulate
   * \param sim_ptrs vector of pointers to the simulation of the fanins.
   * \return A vector of truth-tables, one for each output pin of the node.
   */
  template<typename TT>
  std::vector<TT> compute( node_index_t n, std::vector<TT const*> sim_ptrs ) const
  {
    auto tts = get_functions( n );
    std::vector<TT> res;
    compute( res, n, sim_ptrs );
    return res;
  }

  template<typename TT>
  TT compute( signal_t const& f, std::vector<TT const*> sim_ptrs ) const
  {
    TT res;
    compute( res, f, sim_ptrs );
    return res;
  }

  /*! \brief Inline simulation of the input patterns using the node's function.
   *
   * \param n index of the node to simulate
   * \param sim_ptrs vector of pointers to the simulation of the fanins.
   */
  template<typename TT>
  void compute( std::vector<TT>& res, node_index_t n, std::vector<TT const*> sim_ptrs ) const
  {
    auto simulator_ptr = get_simulator<TT>();
    res.resize( num_outputs( n ) );
    const auto nfanin = fanin_size( n );
    assert( nfanin > 0 );
    assert( sim_ptrs.size() == nfanin );

    _storage->foreach_output_pin( n, [&]( auto const& pin, auto i ) {
      auto id = pin.id;
      auto const& list = _storage->get_chain( id );
      ( *simulator_ptr )( list, sim_ptrs );
      simulator_ptr->get_simulation_inline( res[i], list, sim_ptrs, list.po_at( 0 ) );
    } );
  }

  /*! \brief Inline simulation of the input patterns using the node's function.
   *
   * \param n index of the node to simulate
   * \param sim_ptrs vector of pointers to the simulation of the fanins.
   */
  template<typename TT>
  void compute( TT& res, signal_t const& f, std::vector<TT const*> sim_ptrs ) const
  {
    auto simulator_ptr = get_simulator<TT>();
    const auto nfanin = fanin_size( get_node( f ) );
    assert( nfanin > 0 );
    assert( sim_ptrs.size() == nfanin );
    auto const& g = get_binding( f );
    auto const& list = _storage->get_chain( g.id );
    ( *simulator_ptr )( list, sim_ptrs );
    simulator_ptr->get_simulation_inline( res, list, sim_ptrs, list.po_at( 0 ) );
  }
#pragma endregion

#pragma region Custom node values
  void clear_values() const
  {
    _storage->clear_values();
  }

  uint32_t value( node_index_t n ) const
  {
    return _storage->value( n );
  }

  void set_value( node_index_t n, uint32_t v ) const
  {
    _storage->set_value( n, v );
  }

  uint32_t incr_value( node_index_t n ) const
  {
    return _storage->incr_value( n );
  }

  uint32_t decr_value( node_index_t n ) const
  {
    return _storage->decr_value( n );
  }
#pragma endregion

#pragma region Visited flags
  void clear_visited() const
  {
    _storage->clear_visited();
  }

  auto visited( node_index_t n ) const
  {
    return _storage->visited( n );
  }

  void set_visited( node_index_t n, uint32_t v ) const
  {
    _storage->nodes[n].traversal_id = v;
  }

  uint32_t trav_id() const
  {
    return _storage->get_trav_id();
  }

  void incr_trav_id() const
  {
    _storage->incr_trav_id();
  }

  void set_output_name( uint32_t po_index, std::string const& name )
  {
    _storage->set_output_name( po_index, name );
  }

  void set_input_name( uint32_t pi_index, std::string const& name )
  {
    _storage->set_input_name( pi_index, name );
  }

  bool has_output_name( uint32_t const& po_index ) const
  {
    return _storage->has_output_name( po_index );
  }

  bool has_input_name( uint32_t const& pi_index ) const
  {
    return _storage->has_input_name( pi_index );
  }

  std::string get_output_name( uint32_t const& po_index ) const
  {
    return _storage->get_output_name( po_index );
  }

  std::string get_input_name( uint32_t const& pi_index ) const
  {
    return _storage->get_input_name( pi_index );
  }

  bool has_name( signal_t const& f ) const
  {
    return _storage->has_name( f );
  }

  std::string get_name( signal_t const& f ) const
  {
    return _storage->get_name( f );
  }

  std::string get_network_name() const
  {
    return _storage->get_network_name();
  }

  void set_network_name( std::string const& name )
  {
    _storage->set_network_name( name );
  }

  void set_name( signal_t const& f, std::string const& name )
  {
    _storage->set_name( f, name );
  }
#pragma endregion

#pragma region Getters
  std::vector<signal_t> const& get_children( node_index_t n ) const
  {
    return _storage->get_children( n );
  }

  std::vector<signal_t> const& get_children( signal_t const& f ) const
  {
    return get_children( get_node( f ) );
  }

  std::vector<node_index_t> const& get_fanins( node_index_t n ) const
  {
    return _storage->get_fanins( n );
  }

  double get_area( node_index_t n ) const
  {
    return _storage->get_area( n );
  }

  double get_max_pin_delay( signal_t const& f, uint32_t i ) const
  {
    return _storage->get_max_pin_delay( f, i );
  }

  double get_min_pin_delay( signal_t const& f, uint32_t i ) const
  {
    return _storage->get_min_pin_delay( f, i );
  }

  double get_input_load( signal_t const& f, uint32_t i ) const
  {
    return _storage->get_input_load( f, i );
  }

  auto const& get_library() const
  {
    return _storage->get_library();
  }

  std::vector<kitty::dynamic_truth_table> const get_functions( node_index_t n ) const
  {
    return _storage->get_functions( n );
  }

  std::string const get_expression( signal_t const& f ) const
  {
    auto const n = get_node( f );
    if ( is_pi( n ) )
    {
      std::string expr = "pi" + std::to_string( pi_index( n ) );
      return expr;
    }
    else
      return _storage->get_expression( f );
  }

  uint32_t get_fanin_number( unsigned int id, std::string const& pin_name ) const
  {
    return _storage->get_fanin_number( id, pin_name );
  }
#pragma endregion

#pragma region General methods
  auto& events() const
  {
    return *_events;
  }
#pragma endregion

#pragma region Binding
  std::vector<uint32_t> get_binding_ids( node_index_t n ) const
  {
    return _storage->get_binding_ids( n );
  }

  std::vector<unsigned int> get_binding_ids( std::string const& name ) const
  {
    return _storage->get_binding_ids( name );
  }

  uint32_t get_binding_index( node_index_t n ) const
  {
    return _storage->get_binding_ids( n )[0];
  }

  uint32_t get_binding_index( signal_t const& f ) const
  {
    return _storage->get_binding_ids( f.index )[f.output];
  }

  unsigned int get_pin_id( std::string const& gate_name, std::string const& pin_name ) const
  {
    return _storage->get_pin_id( gate_name, pin_name );
  }

  auto const& get_binding( signal_t const& f ) const
  {
    return _storage->get_binding( f );
  }

  bool has_binding( signal_t const& f ) const
  {
    return _storage->has_binding( f );
  }

  bool has_binding( node_index_t n ) const
  {
    return _storage->has_binding( n );
  }

  bool has_gate( std::string const& name ) const
  {
    return _storage->has_gate( name );
  }

  bool is_input_pin( std::string const& gate_name, std::string const& pin_name ) const
  {
    return _storage->is_input_pin( gate_name, pin_name );
  }

  bool is_output_pin( std::string const& gate_name, std::string const& pin_name ) const
  {
    return _storage->is_output_pin( gate_name, pin_name );
  }
#pragma endregion

public:
  std::shared_ptr<storage_network<DesignType, NumBitsOutputs>> _storage;
  std::shared_ptr<mockturtle::network_events<base_type>> _events;
};

} // namespace network

} // namespace rinox
