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
  \file write_verilog.hpp
  \brief Write networks to structural Verilog format

  \author Andrea Costamagna
*/

#pragma once

#include "../../network/network.hpp"
#include "../../network/signal_map.hpp"
#include "../../traits.hpp"
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/node_map.hpp>
#include <mockturtle/views/topo_view.hpp>

#include <fmt/format.h>
#include <kitty/print.hpp>
#include <lorina/verilog.hpp>

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <regex>
#include <string>
#include <unordered_map>
#include <vector>

namespace rinox
{

namespace io
{

namespace verilog
{
// using namespace std::string_literals;

namespace detail
{

struct bus_info_t
{
  std::string name; // base name, e.g. "A"
  int width;        // number of bits
  bool descending;  // true if bits appear MSB->LSB in the input order
};

class rinox_verilog_writer : public lorina::verilog_writer
{
public:
  rinox_verilog_writer( std::ostream& os )
      : lorina::verilog_writer( os )
  {}

  virtual void on_input( std::vector<bus_info_t> const& inputs ) const
  {
    for ( auto const& bus : inputs )
    {
      if ( bus.width <= 1 )
        _os << fmt::format( "  input {} ;\n", bus.name );
      else
      {
        if ( bus.descending )
          _os << fmt::format( "  input [{}:0] {} ;\n", bus.width - 1, bus.name );
        else
          _os << fmt::format( "  input [0:{}] {} ;\n", bus.width - 1, bus.name );
      }
    }
  }

  virtual void on_output( std::vector<bus_info_t> const& outputs ) const
  {
    for ( auto const& bus : outputs )
    {
      if ( bus.width <= 1 )
        _os << fmt::format( "  output {} ;\n", bus.name );
      else
      {
        if ( bus.descending )
          _os << fmt::format( "  output [{}:0] {} ;\n", bus.width - 1, bus.name );
        else
          _os << fmt::format( "  output [0:{}] {} ;\n", bus.width - 1, bus.name );
      }
    }
  }
};

template<class Ntk>
std::vector<std::pair<bool, std::string>>
format_fanin( Ntk const& ntk, typename Ntk::node const& n, mockturtle::node_map<std::string, Ntk>& node_names )
{
  std::vector<std::pair<bool, std::string>> children;
  ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
    children.emplace_back( std::make_pair( ntk.is_complemented( f ), node_names[f] ) );
  } );
  return children;
}

template<class Ntk>
std::vector<std::pair<bool, std::string>>
format_fanin( Ntk const& ntk, typename Ntk::node const& n, mockturtle::node_map<std::vector<std::string>, Ntk>& node_names )
{
  std::vector<std::pair<bool, std::string>> children;
  ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
    if constexpr ( mockturtle::has_is_multioutput_v<Ntk> )
    {
      children.emplace_back( std::make_pair( ntk.is_complemented( f ), node_names[f][ntk.get_output_pin( f )] ) );
    }
    else
    {
      children.emplace_back( std::make_pair( ntk.is_complemented( f ), node_names[f].front() ) );
    }
  } );
  return children;
}

template<class Ntk>
std::vector<std::pair<bool, std::string>>
format_fanin( Ntk const& ntk, typename Ntk::node const& n, network::incomplete_signal_map<std::string, Ntk>& signal_names )
{
  std::vector<std::pair<bool, std::string>> children;
  ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
    std::string sname = signal_names[f];
    if ( sname.empty() && ntk.is_constant( ntk.get_node( f ) ) )
    {
      sname = signal_names[ntk.get_constant( ntk.constant( ntk.get_node( f ) ) )];
      signal_names[f] = sname;
    }
    children.emplace_back( std::make_pair( ntk.is_complemented( f ), sname ) );
  } );
  return children;
}

template<typename Signal>
struct verilog_writer_signal_hash
{
  uint64_t operator()( const Signal& f ) const
  {
    return f.data;
  }
};

inline std::vector<bus_info_t> infer_buses( const std::vector<std::string>& nets )
{
  // Match "Base[idx]"; if no match, we treat entry as scalar "Base[0]"
  static const std::regex re( R"(([^[]+)\[(\d+)\])" );

  struct Item
  {
    int idx;
    size_t pos;
  }; // bit index and first appearance

  std::unordered_map<std::string, std::vector<Item>> groups;
  groups.reserve( nets.size() );

  for ( size_t pos = 0; pos < nets.size(); ++pos )
  {
    const auto& s = nets[pos];
    std::smatch m;
    std::string base;
    int idx = 0;

    if ( std::regex_match( s, m, re ) )
    {
      base = m[1].str();
      idx = std::stoi( m[2].str() );
    }
    else
    {
      // Plain scalar name -> interpret as base[0]
      base = s;
      idx = 0;
    }

    auto& vec = groups[base];
    // record only the first occurrence position per index
    auto it = std::find_if( vec.begin(), vec.end(), [&]( const Item& it ) { return it.idx == idx; } );
    if ( it == vec.end() )
      vec.push_back( { idx, pos } );
  }

  // Preserve stable order by first appearance of the base name
  auto earliest_pos = []( const std::vector<Item>& items ) {
    size_t p = std::numeric_limits<size_t>::max();
    for ( auto const& it : items )
      p = std::min( p, it.pos );
    return p;
  };

  std::vector<std::pair<std::string, std::vector<Item>>> ordered;
  ordered.reserve( groups.size() );
  for ( auto& kv : groups )
    ordered.emplace_back( kv.first, std::move( kv.second ) );
  std::stable_sort( ordered.begin(), ordered.end(),
                    [&]( auto const& a, auto const& b ) { return earliest_pos( a.second ) < earliest_pos( b.second ); } );

  // Build results
  std::vector<bus_info_t> result;
  result.reserve( ordered.size() );

  for ( auto& [name, items] : ordered )
  {
    if ( items.empty() )
      continue;

    // Sort by first appearance to infer direction from input order
    std::sort( items.begin(), items.end(), []( const Item& a, const Item& b ) { return a.pos < b.pos; } );

    // Extract the seen indices in order of appearance
    std::vector<int> seq;
    seq.reserve( items.size() );
    for ( auto const& it : items )
      seq.push_back( it.idx );

    // Direction: descending if strictly non-increasing and not non-decreasing
    bool nondecreasing = true, nonincreasing = true;
    for ( size_t i = 1; i < seq.size(); ++i )
    {
      if ( seq[i] < seq[i - 1] )
        nondecreasing = false;
      if ( seq[i] > seq[i - 1] )
        nonincreasing = false;
    }
    bool descending = ( !nondecreasing && nonincreasing );

    int width = static_cast<int>( items.size() );

    result.push_back( bus_info_t{ name, width, descending } );
  }

  return result;
}

} // namespace detail

/*! \brief Writes mapped network in structural Verilog format into output stream
 *
 * **Required network functions:**
 * - `num_pis`
 * - `num_pos`
 * - `foreach_pi`
 * - `foreach_node`
 * - `foreach_fanin`
 * - `get_node`
 * - `get_constant`
 * - `is_constant`
 * - `is_pi`
 * - `node_to_index`
 * - `has_binding`
 * - `get_binding_index`
 *
 * \param ntk Mapped network
 * \param os Output stream
 * \param ps Verilog parameters
 */
template<network::design_type_t DesignStyle, uint32_t MaxNumOutputs>
void write_verilog( network::bound_network<DesignStyle, MaxNumOutputs> const& ntk, std::ostream& os, mockturtle::write_verilog_params const& ps = {} )
{
  using Ntk = network::bound_network<DesignStyle, MaxNumOutputs>;
  static_assert( mockturtle::is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( mockturtle::has_num_pis_v<Ntk>, "Ntk does not implement the num_pis method" );
  static_assert( mockturtle::has_num_pos_v<Ntk>, "Ntk does not implement the num_pos method" );
  static_assert( mockturtle::has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( mockturtle::has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( mockturtle::has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
  static_assert( mockturtle::has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( mockturtle::has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( mockturtle::has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( mockturtle::has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( mockturtle::has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( mockturtle::has_has_binding_v<Ntk>, "Ntk does not implement the has_binding method" );
  static_assert( mockturtle::has_get_binding_index_v<Ntk>, "Ntk does not implement the get_binding_index method" );

  assert( ntk.is_combinational() && "Network has to be combinational" );

  rinox::io::verilog::detail::rinox_verilog_writer writer( os );

  std::vector<std::string> inputs;
  if constexpr ( mockturtle::has_has_name_v<Ntk> && mockturtle::has_get_name_v<Ntk> )
  {
    ntk.foreach_pi( [&]( auto const& i, uint32_t index ) {
      auto const f = ntk.make_signal( i );
      if ( ntk.has_name( f ) )
      {
        inputs.emplace_back( ntk.get_name( ntk.make_signal( i ) ) );
      }
      else
      {
        inputs.emplace_back( fmt::format( "x{}", index ) );
      }
    } );
  }
  else
  {
    for ( auto i = 0u; i < ntk.num_pis(); ++i )
    {
      inputs.emplace_back( fmt::format( "x{}", i ) );
    }
  }

  std::vector<std::string> outputs;
  if constexpr ( mockturtle::has_has_output_name_v<Ntk> && mockturtle::has_get_output_name_v<Ntk> )
  {
    ntk.foreach_po( [&]( auto const& o, uint32_t index ) {
      if ( ntk.has_output_name( index ) )
      {
        outputs.emplace_back( ntk.get_output_name( index ) );
      }
      else
      {
        outputs.emplace_back( fmt::format( "y{}", index ) );
      }
    } );
  }
  else
  {
    for ( auto i = 0u; i < ntk.num_pos(); ++i )
    {
      outputs.emplace_back( fmt::format( "y{}", i ) );
    }
  }

  /* compute which nodes are POs and register index */
  network::incomplete_signal_map<std::vector<uint32_t>, Ntk> po_signals( ntk );
  ntk.foreach_po( [&]( auto const& f, auto i ) {
    po_signals[f].push_back( i );
  } );

  std::vector<std::string> ws;
  network::incomplete_signal_map<std::string, Ntk> signal_names( ntk );

  /* constants */
  if ( ntk.has_binding( ntk.get_constant( false ) ) )
  {
    signal_names[ntk.get_constant( false )] = fmt::format( "n{}", ntk.get_constant( false ).index );
    if ( !po_signals.has( ntk.get_constant( false ) ) )
    {
      ws.emplace_back( signal_names[ntk.get_constant( false )] );
    }
  }
  else
  {
    signal_names[ntk.get_constant( false )] = "1'b0";
  }

  if ( ntk.has_binding( ntk.get_constant( true ) ) )
  {
    signal_names[ntk.get_constant( true )] = fmt::format( "n{}", ntk.get_constant( true ).index );
    if ( !po_signals.has( ntk.get_constant( true ) ) )
    {
      ws.emplace_back( signal_names[ntk.get_constant( true )] );
    }
  }
  else
  {
    signal_names[ntk.get_constant( true )] = "1'b1";
  }

  /* add wires */
  ntk.foreach_gate( [&]( auto const& n ) {
    ntk.foreach_output( n, [&]( auto const& f ) {
      if ( !po_signals.has( f ) )
      {
        std::string name = ntk.is_multioutput( n ) ? fmt::format( "n{}_{}", f.index, f.output ) : fmt::format( "n{}", f.index );
        ws.emplace_back( name );
      }
    } );
  } );

  std::string module_name = "top";
  if ( ps.module_name )
  {
    module_name = *ps.module_name;
  }
  else
  {
    if constexpr ( mockturtle::has_get_network_name_v<Ntk> )
    {
      if ( ntk.get_network_name().length() > 0 )
      {
        module_name = ntk.get_network_name();
      }
    }
  }
  auto const info_input = detail::infer_buses( inputs );
  auto const info_output = detail::infer_buses( outputs );
  std::vector<std::string> input_names;
  for ( auto s : info_input )
    input_names.push_back( s.name );
  std::vector<std::string> output_names;
  for ( auto s : info_output )
    output_names.push_back( s.name );

  writer.on_module_begin( module_name, input_names, output_names );

  writer.on_input( info_input );
  writer.on_output( info_output );

  if ( !ws.empty() )
  {
    writer.on_wire( ws );
  }

  ntk.foreach_pi( [&]( auto const& n, auto i ) {
    signal_names[ntk.make_signal( n )] = inputs[i];
  } );

  auto const& gates = ntk.get_library();

  int nDigits = (int)std::floor( std::log10( ntk.num_gates() + 2 + ntk.num_pos() ) );
  unsigned int length = 0;
  unsigned counter = 0;

  for ( auto const& gate : gates )
  {
    length = std::max( length, static_cast<unsigned int>( gate.name.length() ) );
  }

  mockturtle::topo_view ntk_topo{ ntk };
  std::vector<std::pair<std::string, std::string>> assignments;

  ntk_topo.foreach_node( [&]( auto const& n ) {
    if ( ntk_topo.is_constant( n ) && ( n > 1 ) )
      return true;
    ntk_topo.foreach_output( n, [&]( auto const& f ) {
      if ( po_signals.has( f ) && !ntk_topo.is_constant( ntk_topo.get_node( f ) ) )
      {
        if ( ntk.has_name( f ) )
        {
          auto const sname = ntk.get_name( f );
          if ( !sname.empty() )
            signal_names[f] = sname;
          else
            signal_names[f] = outputs[po_signals[f][0]];
          auto const nout = po_signals[f].size();
          for ( auto i = 0u; i < nout; ++i )
          {
            auto const oname = ntk.get_output_name( po_signals[f][i] );
            if ( oname != signal_names[f] )
            {
              assignments.emplace_back( std::make_pair( oname, signal_names[f] ) );
            }
          }
        }
        else
          signal_names[f] = outputs[po_signals[f][0]];
      }
      else if ( !ntk.is_constant( n ) && !ntk.is_pi( n ) )
      {
        if ( ntk.has_name( f ) )
          signal_names[f] = ntk.get_name( f );
        else
          signal_names[f] = ntk.is_multioutput( n ) ? fmt::format( "n{}_{}", f.index, f.output ) : fmt::format( "n{}", f.index );
      }
    } );

    if ( ntk.has_binding( n ) )
    {
      auto const& gate = gates[ntk.get_binding_index( n )];
      std::string name = gate.name;

      int digits = counter == 0 ? 0 : (int)std::floor( std::log10( counter ) );
      auto fanin_names = detail::format_fanin<Ntk>( ntk, n, signal_names );
      std::vector<std::pair<std::string, std::string>> args;

      auto i = 0;
      for ( auto pair : fanin_names )
      {
        args.emplace_back( std::make_pair( gate.pins[i++].name, pair.second ) );
      }
      ntk.foreach_output( n, [&]( auto const& f ) {
        args.emplace_back( std::make_pair( gates[ntk.get_binding_index( f )].output_name, signal_names[f] ) );
      } );

      writer.on_module_instantiation( name.append( std::string( length - name.length(), ' ' ) ),
                                      {},
                                      std::string( "g" ) + std::string( nDigits - digits, '0' ) + std::to_string( counter ),
                                      args );
      ++counter;

      /* if node drives multiple POs, duplicate */
      ntk.foreach_output( n, [&]( auto const& f ) {
        if ( po_signals.has( f ) && po_signals[f].size() > 1 )
        {
          if ( ps.verbose )
          {
            std::cerr << "[i] signal {" << f.index << ", " << f.output << "} driving multiple POs has been duplicated.\n";
          }
          auto const& po_list = po_signals[f];
          for ( auto i = 1u; i < po_list.size(); ++i )
          {
            digits = counter == 0 ? 0 : (int)std::floor( std::log10( counter ) );
            args[args.size() - 1] = std::make_pair( gates[ntk.get_binding_index( f )].output_name, outputs[po_list[i]] );

            writer.on_module_instantiation( name.append( std::string( length - name.length(), ' ' ) ),
                                            {},
                                            std::string( "g" ) + std::string( nDigits - digits, '0' ) + std::to_string( counter ),
                                            args );
            ++counter;
          }
        }
      } );
    }
    else if ( !ntk.is_constant( n ) && !ntk.is_pi( n ) )
    {
      std::cerr << "[e] internal node " << n << " is not mapped.\n";
    }

    return true;
  } );

  kitty::dynamic_truth_table inv( 1u );
  kitty::create_from_binary_string( inv, "01" );
  auto id = ntk_topo.get_id( inv );
  auto const& gate_inv = gates[id];

  ntk_topo.foreach_po( [&]( auto const& f, auto io ) {
    auto const n = ntk_topo.get_node( f );
    if ( ntk_topo.is_constant( n ) )
    {
      int digits = counter == 0 ? 0 : (int)std::floor( std::log10( counter ) );
      std::vector<std::pair<std::string, std::string>> args;
      if ( ntk_topo.constant( n ) )
      {
        args.emplace_back( std::make_pair( gate_inv.output_name, outputs[io] ) );
        args.emplace_back( std::make_pair( gate_inv.pins[0].name, signal_names[ntk_topo.get_constant( false )] ) );
      }
      else
      {
        args.emplace_back( std::make_pair( gate_inv.output_name, outputs[io] ) );
        args.emplace_back( std::make_pair( gate_inv.pins[0].name, signal_names[ntk_topo.get_constant( true )] ) );
      }
      std::string name = gate_inv.name;

      writer.on_module_instantiation( name.append( std::string( length - name.length(), ' ' ) ),
                                      {},
                                      std::string( "g" ) + std::string( nDigits - digits, '0' ) + std::to_string( counter ),
                                      args );
      ++counter;
    }
  } );

  std::string lhs;
  std::vector<std::pair<bool, std::string>> rhs;
  std::string const op = "";
  for ( auto a : assignments )
  {
    lhs = a.first;
    rhs = { std::make_pair( false, a.second ) };
    writer.on_assign( lhs, rhs, op );
  }

  writer.on_module_end();
}

/*! \brief Writes network in structural Verilog format into a file
 *
 * \param ntk Network
 * \param filename Filename
 */
template<network::design_type_t DesignStyle, uint32_t MaxNumOutputs>
void write_verilog( network::bound_network<DesignStyle, MaxNumOutputs> const& ntk, std::string const& filename, mockturtle::write_verilog_params const& ps = {} )
{
  std::ofstream os( filename.c_str(), std::ofstream::out );
  write_verilog<DesignStyle, MaxNumOutputs>( ntk, os, ps );
  os.close();
}

} /* namespace verilog */

} /* namespace io */

} /* namespace rinox */