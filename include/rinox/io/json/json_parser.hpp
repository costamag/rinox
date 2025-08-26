/* rinox: C++ parsing library
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
  \file json_parser.hpp
  \brief Implements simplistic Json parser

  \author Andrea Costamagna
*/

#pragma once
#include <fmt/format.h>
#include <optional>
#include <queue>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <lorina/detail/tokenizer.hpp>
#include <lorina/detail/utils.hpp>
#include <lorina/diagnostics.hpp>

#include <rinox/diagnostics/diagnostics.hpp>
#include "rinox/io/json/json_stream.hpp"

namespace rinox
{

namespace io
{

namespace json
{

template<typename reader_t>
class json_parser
{
public:
  using ModuleInstFn = lorina::detail::Func<
      std::string,                                     // module_name
      std::vector<std::string>,                        // params (if any)
      std::string,                                     // instance_name
      std::vector<std::pair<std::string, std::string>> // args: ".pin" -> net
      >;

  using CellFn = lorina::detail::Func<
      std::vector<std::pair<std::string, std::string>>, // input map
      std::vector<std::pair<std::string, std::string>>, // output map
      std::vector<unsigned int>                         // binding ids
      >;

  using ModuleInstParamMap = lorina::detail::ParamPackMap<
      std::string, // key
      std::string,
      std::vector<std::string>,
      std::string,
      std::vector<std::pair<std::string, std::string>>>;

  using CellParamMap = lorina::detail::ParamPackMap<
      std::string,                                      // key
      std::vector<std::pair<std::string, std::string>>, // input map
      std::vector<std::pair<std::string, std::string>>, // output map
      std::vector<unsigned int>>;

  using PackedFns = lorina::detail::FuncPackN<ModuleInstFn, CellFn>;
  using ParamMaps = lorina::detail::ParamPackMapN<ModuleInstParamMap, CellParamMap>;

  static constexpr int MODULE_INST_FN = 0;
  static constexpr int CELL_FN = 1;

public:
  json_parser( std::istream& in,
               const reader_t& rdr,
               lorina::diagnostic_engine* de = nullptr,
               std::string_view module = "top" )
      : jstream_( in, de ), reader_( rdr ), diag_( de ), on_action_( PackedFns(
                                                             ModuleInstFn( [&]( const std::string mn,
                                                                                const std::vector<std::string>& params,
                                                                                const std::string inst,
                                                                                const std::vector<std::pair<std::string, std::string>>& pin2pin ) {
                                                               reader_.on_module_instantiation( mn, params, inst, pin2pin );
                                                             } ),
                                                             CellFn( [&]( const std::vector<std::pair<std::string, std::string>>& in_map,
                                                                          const std::vector<std::pair<std::string, std::string>>& out_map,
                                                                          const std::vector<unsigned int>& ids ) {
                                                               reader_.on_cell( in_map, out_map, ids );
                                                             } ) ) )
  {
    on_action_.declare_known( "0" );
    on_action_.declare_known( "1" );
    on_action_.declare_known( "x" );
    on_action_.declare_known( "z" );

    for ( auto& n : jstream_.module_names() )
      modules_.insert( n );
  }

  bool parse_modules()
  {
    for ( auto const& name : jstream_.module_names() )
    {
      if ( !parse_module( name ) )
      {
        rinox::diagnostics::REPORT_DIAG( diag_, lorina::diagnostic_level::fatal,
                     "failed to parse module `{}`", name.c_str() );
        return false;
      }
    }
    return true;
  }

private:
  static std::string size_from_bits_( size_t n )
  {
    return ( n <= 1 ) ? std::string{} : fmt::format( "{}:0", n - 1 );
  }

  static std::string bit_token_( const rinox::io::json::json_bit_t& b )
  {
    if ( std::holds_alternative<int64_t>( b.v ) )
      return std::to_string( std::get<int64_t>( b.v ) );
    return std::get<std::string>( b.v ); // "0","1","x","z"
  }

  bool parse_module( std::string const& name )
  {
    if ( !jstream_.set_module( name ) )
    {
      rinox::diagnostics::REPORT_DIAG( diag_, lorina::diagnostic_level::fatal,
                   "module `{}` not found in JSON", name.c_str() );
      return false;
    }

    std::vector<rinox::io::json::port_instance_t> ports;
    std::queue<rinox::io::json::instance_t> pending; // cells/netnames buffered

    rinox::io::json::instance_t inst;
    for ( ;; )
    {
      auto rc = jstream_.get_instance( inst );
      if ( rc == rinox::io::json::instance_return_code::invalid )
      {
        rinox::diagnostics::REPORT_DIAG( diag_, lorina::diagnostic_level::fatal,
                     "failed to parse a new instance" );
        return false;
      }
      if ( rc == rinox::io::json::instance_return_code::end )
        break;

      if ( std::holds_alternative<rinox::io::json::port_instance_t>( inst ) )
        ports.push_back( std::get<rinox::io::json::port_instance_t>( inst ) );
      else
        pending.push( std::move( inst ) );
    }

    std::vector<std::string> inouts;
    inouts.reserve( ports.size() );
    for ( auto const& p : ports )
      inouts.push_back( p.name );
    reader_.on_module_header( name, inouts );

    for ( auto const& p : ports )
    {
      const auto sz = size_from_bits_( p.bits.size() );
      if ( p.direction == "input" )
      {
        std::vector<std::string> name_ids;
        for ( auto const& b : p.bits )
        {
          if ( std::holds_alternative<int64_t>( b.v ) )
          {
            const std::string bit_id = std::to_string( std::get<int64_t>( b.v ) );
            name_ids.push_back( bit_id );
            on_action_.declare_known( bit_id );
          }
          else
            on_action_.declare_known( std::get<std::string>( b.v ) ); // "0","1","x","z"
        }
        reader_.on_input( p.name, name_ids, p.upto );
      }
      else if ( p.direction == "output" )
      {
        std::vector<std::string> name_ids;
        for ( auto const& b : p.bits )
        {
          if ( std::holds_alternative<int64_t>( b.v ) )
          {
            const std::string bit_id = std::to_string( std::get<int64_t>( b.v ) );
            name_ids.push_back( bit_id );
          }
          else
          {
            auto const s = std::get<std::string>( b.v );
            name_ids.push_back( s );
          }
        }
        reader_.on_output( p.name, name_ids, p.upto );
      }
      else
      {
        rinox::diagnostics::REPORT_DIAG( diag_, lorina::diagnostic_level::fatal,
                     "port direction `{}` not supported", p.direction.c_str() );
        return false;
      }
    }

    while ( !pending.empty() )
    {
      auto cur = std::move( pending.front() );
      pending.pop();

      if ( std::holds_alternative<rinox::io::json::cell_instance_t>( cur ) )
      {
        const auto& cell = std::get<rinox::io::json::cell_instance_t>( cur );

        std::vector<std::pair<std::string, std::string>> inmap, outmap;
        std::vector<std::string> dep_inputs, dep_outputs;

        auto is_output_pin = [&]( std::string const& pin ) -> bool {
          auto it = cell.port_dirs.find( pin );
          if ( it != cell.port_dirs.end() )
            return it->second == "output";
          return reader_.is_output_pin( cell.type, pin );
        };

        for ( auto const& kv : cell.connections )
        {
          const auto& pin = kv.first;
          const auto& bits = kv.second;
          if ( bits.empty() )
            continue;

          const std::string net = bit_token_( bits.front() );

          if ( is_output_pin( pin ) )
          {
            outmap.emplace_back( pin, net );
            dep_outputs.push_back( net );
          }
          else
          {
            inmap.emplace_back( pin, net );
            dep_inputs.push_back( net );
          }
        }

        if ( modules_.count( cell.type ) )
        {
          std::vector<std::pair<std::string, std::string>> args;
          args.reserve( inmap.size() + outmap.size() );
          for ( auto& kv : inmap )
            args.emplace_back( "." + kv.first, kv.second );
          for ( auto& kv : outmap )
            args.emplace_back( "." + kv.first, kv.second );

          std::vector<std::string> params;

          on_action_.call_deferred<MODULE_INST_FN>(
              /* deps in  */ dep_inputs,
              /* deps out */ dep_outputs,
              /* params   */ std::make_tuple( cell.type, params, cell.name, args ) );
        }
        else
        {
          // Primitive / tech cell
          std::vector<unsigned int> ids;
          for ( auto dep : outmap )
            ids.push_back( reader_.get_pin_id( cell.type, dep.first ) );

          on_action_.call_deferred<CELL_FN>(
              /* deps in  */ dep_inputs,
              /* deps out */ dep_outputs,
              /* params   */ std::make_tuple( inmap, outmap, ids ) );
        }
      }
    }

    bool ok = true;
    const auto& deps = on_action_.unresolved_dependencies();
    if ( !deps.empty() )
      ok = false;

    for ( const auto& r : deps )
    {
      if ( diag_ )
      {
        diag_->report( lorina::diag_id::WRN_UNRESOLVED_DEPENDENCY )
            .add_argument( r.first )
            .add_argument( r.second );
      }
    }

    if ( !ok )
    {
      rinox::diagnostics::REPORT_DIAG( diag_, lorina::diagnostic_level::note,
                   "dangling objects not parsed" );
      return false;
    }

    reader_.on_endmodule();
    return true;
  }

private:
  rinox::io::json::json_stream jstream_;
  const reader_t& reader_;
  lorina::diagnostic_engine* diag_;
  lorina::detail::call_in_topological_order<PackedFns, ParamMaps> on_action_;
  std::unordered_set<std::string> modules_;
};

} // namespace json

} // namespace io

} // namespace rinox
