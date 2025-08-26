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
  \file augmented_library.hpp
  \brief Implements methods for handling and evaluating a library of standard cells.

  This engine can be used for efficient Boolean evaluation of the gates in a
  standard cell library. Each gate is represented as an AIG index chain for
  efficient evaluation. Additionally, in the presence of multiple-output cells,
  this engine identifies which gates belong to a multiple-output and allows handling
  this information.

  NOTE: The augmented library can be made arbitrarily complex adding technological
  information. This data structure can be modified to store detailed information from
  the liberty file.

  \author Andrea Costamagna
*/

#pragma once

#include <set>
#include <vector>

#include "../synthesis/xaig_decompose.hpp"

#include <kitty/print.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/truth_table_cache.hpp>

namespace rinox
{

namespace libraries
{

template<network::design_type_t DesignType>
class augmented_library;

template<>
class augmented_library<network::design_type_t::CELL_BASED>
{
public:
  /* truth table used to express the gate functionality */
  using func_t = kitty::dynamic_truth_table;
  using chain_t = evaluation::chains::large_xag_chain;
  /* gate type in the raw format */
  using gate = mockturtle::gate;

  /*! \brief Augmented gate.
   *
   * A raw gate is augmented by decomposing it into an index chain for
   * efficient simulation.
   */
  struct gate_t : gate
  {
    chain_t aig_chain;
    std::vector<double> max_pin_time;
    std::vector<double> min_pin_time;
    double avg_pin_delay;

    gate_t( const gate& g, chain_t const& chain )
        : gate( g ),
          aig_chain( chain ),
          max_pin_time( g.num_vars, std::numeric_limits<double>::min() ),
          min_pin_time( g.num_vars, std::numeric_limits<double>::max() )
    {
      avg_pin_delay = 0;
      for ( auto i = 0u; i < g.num_vars; ++i )
      {
        double const rise_time = g.pins[i].rise_block_delay;
        double const fall_time = g.pins[i].fall_block_delay;
        max_pin_time[i] = std::max( rise_time, fall_time );
        min_pin_time[i] = std::min( rise_time, fall_time );
        avg_pin_delay += 0.5 * ( max_pin_time[i] + min_pin_time[i] );
      }
      avg_pin_delay /= (double)g.num_vars;
    }
  };

  /*! \brief Construction via specification of the simpler library.
   *
   * The gates should specify at least the gate's functionality, from
   * which this constructor can synthesize an index chain for each gate.
   */
  augmented_library( std::vector<gate> const& raw_gates )
      : synth( st ),
        raw_gates_( raw_gates )
  {
    gates_.reserve( raw_gates.size() );
    for ( gate const& g : raw_gates )
    {
      add_gate( g );
      name_to_ids_[g.name].push_back( g.id );
    }

    for ( gate const& g : raw_gates )
    {
      auto it = single_output.find( g.name );
      if ( it != single_output.end() )
      {
        tt_to_index_[g.function] = g.id;
      }
    }
  }

  bool constant( std::vector<unsigned int> const& ids ) const
  {
    bool cst = false;
    for ( auto const& id : ids )
      cst |= ( kitty::count_ones(raw_gates_[id].function) > 0 );
    return cst;
  }

  bool is_constant( std::vector<unsigned int> const& ids ) const
  {
    bool constant = false;
    for ( auto const& id : ids )
      constant |= (raw_gates_[id].function.num_vars() == 0);
    return constant;
  }

  /*! \brief Augment the gate and add it to the library.
   *
   * \param raw_gate Representation containing the gate's function
   */
  uint32_t add_gate( gate const& g )
  {
    synth( g.function );
    chain_t const chain = synth.get_chain();

    /* check if the node should be added to the multiple output gates */
    if ( single_output.find( g.name ) != single_output.end() )
    {
      single_output.erase( g.name );
      multiple_output.insert( g.name );
    }
    else
    {
      single_output.insert( g.name );
    }

    auto binding_id = static_cast<uint32_t>( gates_.size() );
    gates_.emplace_back( g, chain );
    return binding_id;
  }

  /*! \brief Getter of the chain synthesizing the gate's functionality. */
  chain_t const& get_chain( uint32_t id ) const
  {
    return gates_[id].aig_chain;
  }

  /*! \brief Getter of the gate's name. */
  std::string const& get_name( uint32_t id ) const
  {
    return gates_[id].name;
  }

  /*! \brief Getter of the augmented gate. */
  gate_t const& get_gate( uint32_t id ) const
  {
    return gates_[id];
  }

  /*! \brief Getter of the gate's area. */
  double const& get_area( uint32_t id ) const
  {
    return gates_[id].area;
  }

  double get_max_pin_delay( uint32_t id, uint32_t i ) const
  {
    return gates_[id].max_pin_time[i];
  }

  double get_min_pin_delay( uint32_t id, uint32_t i ) const
  {
    return gates_[id].min_pin_time[i];
  }

  double get_input_load( uint32_t id, uint32_t i ) const
  {
    return gates_[id].pins[i].input_load;
  }

  std::vector<gate> const& get_raw_gates() const
  {
    return raw_gates_;
  }

  std::vector<gate_t> const& get_aug_gates() const
  {
    return gates_;
  }

  std::vector<unsigned int> get_binding_ids( std::string const& gate_name ) const
  {
    auto const it = name_to_ids_.find( gate_name );
    if ( it == name_to_ids_.end() )
    {
      std::cerr << "[e] No binding found for the gate " << gate_name << std::endl;
      return {};
    }
    return it->second;
  }

  unsigned int get_pin_id( std::string const& gate_name, std::string const& pin_name ) const
  {
    auto const ids = get_binding_ids( gate_name );
    for ( auto id : ids )
    {
      if ( raw_gates_[id].output_name == pin_name )
        return id;
    }
    std::cerr << "[e] pin name " << pin_name << " not found in " << gate_name << std::endl;
    return std::numeric_limits<unsigned int>::max();
  }

  std::optional<unsigned int> get_id( kitty::dynamic_truth_table const& tt ) const
  {
    auto it = tt_to_index_.find( tt );
    if ( it != tt_to_index_.end() )
    {
      return it->second;
    }
    return std::nullopt;
  }

  uint32_t get_fanin_number( unsigned int id, std::string const& pin_name ) const
  {
    auto const& g = gates_[id];
    for ( uint32_t i = 0; i < g.num_vars; ++i )
    {
      if ( g.pins[i].name == pin_name )
      {
        return i;
      }
    }
    std::cerr << "[e] Pin " << pin_name << " not found in gate " << g.name << std::endl;
    return std::numeric_limits<uint32_t>::max();
  }

  /*! \brief Check if the gate is a multiple output gate from its name */
  bool is_multioutput( std::string const& name ) const
  {
    return multiple_output.find( name ) != multiple_output.end();
  }

  bool has_gate( std::string const& name ) const
  {
    return single_output.find( name ) != single_output.end() ||
           multiple_output.find( name ) != multiple_output.end();
  }

  bool is_input_pin( std::string const& gate_name, std::string const& pin_name ) const
  {
    auto const it = name_to_ids_.find( gate_name );
    if ( it == name_to_ids_.end() )
    {
      return false;
    }
    for ( auto const& id : it->second )
    {
      auto const& g = gates_[id];
      for ( auto const& pin : g.pins )
      {
        if ( pin.name == pin_name )
        {
          return true;
        }
      }
    }
    return false;
  }

  bool is_output_pin( std::string const& gate_name, std::string const& pin_name ) const
  {
    auto const it = name_to_ids_.find( gate_name );
    if ( it == name_to_ids_.end() )
    {
      return false;
    }
    for ( auto const& id : it->second )
    {
      auto const& g = gates_[id];
      if ( g.output_name == pin_name )
        return true;
    }
    return false;
  }

private:
  /*! \brief Augmented technology library */
  std::vector<gate> raw_gates_;
  std::vector<gate_t> gates_;
  /*! \brief Synthesis engine for AIG index chains */
  synthesis::xaig_decompose_stats st;
  synthesis::xaig_decompose<false, false> synth;
  /* contains the name of the multiple-output gates in the library */
  std::set<std::string> multiple_output;
  std::set<std::string> single_output;
  phmap::flat_hash_map<kitty::dynamic_truth_table, unsigned int, kitty::hash<kitty::dynamic_truth_table>> tt_to_index_;
  std::unordered_map<std::string, std::vector<unsigned int>> name_to_ids_;
};

template<>
class augmented_library<network::design_type_t::ARRAY_BASED>
{
public:
  /* truth table used to express the gate functionality */
  using func_t = kitty::dynamic_truth_table;
  using chain_t = evaluation::chains::large_xag_chain;
  /* gate type in the raw format */
  struct gate
  {
    gate( kitty::dynamic_truth_table const& function )
        : function( function )
    {}

    kitty::dynamic_truth_table function;
  };

  struct gate_t : gate
  {
    gate_t( kitty::dynamic_truth_table const& function, chain_t const& chain )
        : gate( function ), aig_chain( chain )
    {}

    chain_t aig_chain;
  };

  /*! \brief Construction via specification of the simpler library.
   *
   * The gates should specify at least the gate's functionality, from
   * which this constructor can synthesize an index chain for each gate.
   */
  augmented_library( uint32_t capacity = 1000 )
      : synth( st )
  {
    tt_to_index_.reserve( capacity );
    gates_.reserve( capacity );
  }

  /*! \brief Augment the gate and add it to the library.
   *
   * \param raw_gate Representation containing the gate's function
   */
  uint32_t add_gate( kitty::dynamic_truth_table const& function )
  {
    /* is truth table already in cache? */
    const auto it = tt_to_index_.find( function );
    if ( it != tt_to_index_.end() )
    {
      return it->second;
    }
    synth( function );
    chain_t const chain = synth.get_chain();
    auto const binding_id = static_cast<uint32_t>( gates_.size() );
    tt_to_index_[function] = binding_id;
    gates_.emplace_back( gate_t{ function, chain } );
    return binding_id;
  }

  /*! \brief Augment the gate and add it to the library.
   *
   * \param raw_gate Representation containing the gate's function
   */
  std::vector<uint32_t> add_gates( std::vector<kitty::dynamic_truth_table> const& functions )
  {
    std::vector<uint32_t> binding_ids;
    for ( auto const& func : functions )
    {
      binding_ids.push_back( add_gate( func ) );
    }
    return binding_ids;
  }

  /*! \brief Getter of the chain synthesizing the gate's functionality. */
  chain_t const& get_chain( uint32_t id ) const
  {
    return gates_[id].aig_chain;
  }

  /*! \brief Getter of the gate's name. */
  std::string const get_name( uint32_t id ) const
  {
    return kitty::to_hex( gates_[id].function );
  }

  /*! \brief Getter of the gate's area. */
  double const get_area( uint32_t id ) const
  {
    (void)id;
    return 1.0;
  }

  /*! \brief Getter of the augmented gate. */
  gate_t const& get_gate( uint32_t id ) const
  {
    return gates_[id];
  }

  double get_max_pin_delay( uint32_t id, uint32_t i ) const
  {
    (void)id;
    (void)i;
    return 1.0;
  }

  double get_min_pin_delay( uint32_t id, uint32_t i ) const
  {
    (void)id;
    (void)i;
    return 1.0;
  }

  double get_input_load( uint32_t id, uint32_t i ) const
  {
    (void)id;
    (void)i;
    return 1.0;
  }

  std::vector<gate_t> const& get_aug_gates() const
  {
    return gates_;
  }

  std::optional<unsigned int> get_id( kitty::dynamic_truth_table const& tt ) const
  {
    auto it = tt_to_index_.find( tt );
    if ( it != tt_to_index_.end() )
    {
      return it->second;
    }
    return std::nullopt;
  }

  /*! \brief Check if the gate is a multiple output gate from its name */
  bool is_multioutput( std::string const& name ) const
  {
    (void)name;
    return false;
  }

  bool constant( std::vector<unsigned int> const& ids ) const
  {
    bool cst = false;
    for ( auto const& id : ids )
      cst |= ( kitty::count_ones(gates_[id].function) > 0 );
    return cst;
  }

  bool is_constant( std::vector<unsigned int> const& ids ) const
  {
    bool constant;
    for ( auto const& id : ids )
      constant |= gates_[id].function.num_vars() == 0 ;
    return constant;
  }

private:
  /*! \brief Augmented technology library */
  std::vector<gate_t> gates_;
  phmap::flat_hash_map<kitty::dynamic_truth_table, uint32_t, kitty::hash<kitty::dynamic_truth_table>> tt_to_index_;

  /*! \brief Synthesis engine for AIG index chains */
  synthesis::xaig_decompose_stats st;
  synthesis::xaig_decompose<false, false> synth;
};

} // namespace libraries

} // namespace rinox