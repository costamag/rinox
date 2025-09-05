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
  \file json_stream.hpp
  \brief Core engine for FSM-like parsing of json instances

  \author Andrea Costamagna
*/

#pragma once
#include <rinox/diagnostics/diagnostics.hpp>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <istream>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

namespace rinox
{

namespace io
{

namespace json
{

struct json_bit_t
{
  std::variant<int64_t, std::string> v;
};

inline bool operator==( const json_bit_t& a, const json_bit_t& b ) noexcept
{
  return a.v == b.v;
}
inline bool operator!=( const json_bit_t& a, const json_bit_t& b ) noexcept
{
  return !( a == b );
}

struct port_instance_t
{
  std::string name, direction;
  std::vector<json_bit_t> bits;
  int64_t offset = 0;
  bool upto = false;
  bool is_signed = false;
  int hide_name = 0;
  std::vector<std::pair<std::string, std::string>> attributes;
};

struct cell_instance_t
{
  std::string name, type;
  std::string model;
  std::unordered_map<std::string, std::string> port_dirs;
  std::unordered_map<std::string, std::vector<json_bit_t>> connections;
  std::vector<std::pair<std::string, std::string>> parameters;
  std::vector<std::pair<std::string, std::string>> attributes;
};

struct net_name_instance_t
{
  std::string name;
  std::vector<json_bit_t> bits;
};

using instance_t = std::variant<port_instance_t, cell_instance_t, net_name_instance_t>;
enum class instance_return_code
{
  valid,
  invalid,
  end
};

class json_stream
{
public:
  explicit json_stream( std::istream& in, lorina::diagnostic_engine* diag = nullptr )
      : diag_( diag )
  {
    rapidjson::IStreamWrapper isw( in );
    doc_.ParseStream( isw );
    if ( doc_.HasParseError() || !doc_.IsObject() )
    {
      parse_ok_ = false;
      return;
    }
    if ( !doc_.HasMember( "modules" ) || !doc_["modules"].IsObject() )
    {
      parse_ok_ = false;
      return;
    }

    std::vector<std::string> mnames = module_names();
    assert( ( mnames.size() == 1 ) && "only stable for one module" );
    module_name_ = mnames[0];
    set_module( module_name_ );
  }

  std::vector<std::string> module_names() const
  {
    std::vector<std::string> names;
    if ( !doc_.IsObject() || !doc_.HasMember( "modules" ) || !doc_["modules"].IsObject() )
      return names;
    const auto& mods = doc_["modules"];
    for ( auto it = mods.MemberBegin(); it != mods.MemberEnd(); ++it )
      names.push_back( it->name.GetString() );
    return names;
  }

  bool set_module( std::string_view name )
  {
    if ( !doc_.IsObject() || !doc_.HasMember( "modules" ) || !doc_["modules"].IsObject() )
      return false;
    const auto& mods = doc_["modules"];
    for ( auto it = mods.MemberBegin(); it != mods.MemberEnd(); ++it )
    {
      if ( name == it->name.GetString() )
      {
        module_ = &it->value;
        module_name_ = it->name.GetString();

        ports_ = cells_ = nets_ = nullptr;
        if ( module_->HasMember( "ports" ) && ( *module_ )["ports"].IsObject() )
        {
          ports_ = &( *module_ )["ports"];
          p_it_ = ports_->MemberBegin();
          p_end_ = ports_->MemberEnd();
        }
        if ( module_->HasMember( "cells" ) && ( *module_ )["cells"].IsObject() )
        {
          cells_ = &( *module_ )["cells"];
          c_it_ = cells_->MemberBegin();
          c_end_ = cells_->MemberEnd();
        }
        if ( module_->HasMember( "netnames" ) && ( *module_ )["netnames"].IsObject() )
        {
          nets_ = &( *module_ )["netnames"];
          n_it_ = nets_->MemberBegin();
          n_end_ = nets_->MemberEnd();
        }
        return true;
      }
    }
    return false;
  }

  instance_return_code get_instance( instance_t& out )
  {
    if ( !parse_ok_ )
      return instance_return_code::invalid;

    try
    {
      if ( ports_ && p_it_ != p_end_ )
        return read_port_( out );
      if ( cells_ && c_it_ != c_end_ )
        return read_cell_( out );
      if ( nets_ && n_it_ != n_end_ )
        return read_net_( out );
      return instance_return_code::end;
    }
    catch ( ... )
    {
      return instance_return_code::invalid;
    }
  }

  const std::string& module_name() const { return module_name_; }
  int file_line = 0;

private:
  static bool is_number_( const rapidjson::Value& v )
  {
    return v.IsInt64() || v.IsUint64() || v.IsInt() || v.IsUint();
  }
  static int64_t to_int64_( const rapidjson::Value& v )
  {
    if ( v.IsInt64() )
      return v.GetInt64();
    if ( v.IsUint64() )
      return static_cast<int64_t>( v.GetUint64() );
    if ( v.IsInt() )
      return v.GetInt();
    if ( v.IsUint() )
      return static_cast<int64_t>( v.GetUint() );
    return 0;
  }

  static std::vector<json_bit_t> read_bits_( const rapidjson::Value& v )
  {
    std::vector<json_bit_t> out;
    if ( !v.IsArray() )
      return out;
    for ( auto it = v.Begin(); it != v.End(); ++it )
    {
      if ( is_number_( *it ) )
      {
        out.push_back( json_bit_t{ to_int64_( *it ) } );
      }
      else if ( it->IsString() )
      {
        out.push_back( json_bit_t{ std::string( it->GetString(), it->GetStringLength() ) } );
      }
    }
    return out;
  }

  static std::string raw_to_string_( const rapidjson::Value& v )
  {
    rapidjson::StringBuffer sb;
    rapidjson::Writer<rapidjson::StringBuffer> w( sb );
    v.Accept( w );
    return std::string( sb.GetString(), sb.GetSize() );
  }

  instance_return_code read_port_( instance_t& out )
  {
    const auto& mem = *p_it_;
    const std::string name = mem.name.GetString();
    const auto& pobj = mem.value;
    if ( !pobj.IsObject() )
    {
      ++p_it_;
      return instance_return_code::invalid;
    }

    port_instance_t port;
    port.name = name;

    if ( pobj.HasMember( "direction" ) && pobj["direction"].IsString() )
      port.direction = pobj["direction"].GetString();
    if ( pobj.HasMember( "bits" ) )
      port.bits = read_bits_( pobj["bits"] );

    if ( pobj.HasMember( "offset" ) && is_number_( pobj["offset"] ) )
      port.offset = to_int64_( pobj["offset"] );
    if ( pobj.HasMember( "upto" ) )
      port.upto = pobj["upto"].IsBool() ? pobj["upto"].GetBool() : ( is_number_( pobj["upto"] ) && to_int64_( pobj["upto"] ) != 0 );
    if ( pobj.HasMember( "signed" ) )
      port.is_signed = pobj["signed"].IsBool() ? pobj["signed"].GetBool() : ( is_number_( pobj["signed"] ) && to_int64_( pobj["signed"] ) != 0 );
    if ( pobj.HasMember( "hide_name" ) && is_number_( pobj["hide_name"] ) )
      port.hide_name = static_cast<int>( to_int64_( pobj["hide_name"] ) );
    if ( pobj.HasMember( "attributes" ) && pobj["attributes"].IsObject() )
    {
      for ( auto it = pobj["attributes"].MemberBegin(); it != pobj["attributes"].MemberEnd(); ++it )
      {
        port.attributes.emplace_back( it->name.GetString(), raw_to_string_( it->value ) );
      }
    }

    ++p_it_;
    out = std::move( port );
    return instance_return_code::valid;
  }

  instance_return_code read_cell_( instance_t& out )
  {
    const auto& mem = *c_it_;
    const std::string name = mem.name.GetString();
    const auto& cobj = mem.value;
    if ( !cobj.IsObject() )
    {
      ++c_it_;
      return instance_return_code::invalid;
    }

    cell_instance_t cell;
    cell.name = name;
    if ( cobj.HasMember( "type" ) && cobj["type"].IsString() )
      cell.type = cobj["type"].GetString();

    if ( cobj.HasMember( "model" ) && cobj["model"].IsString() )
      cell.model = cobj["model"].GetString();

    if ( cobj.HasMember( "port_directions" ) && cobj["port_directions"].IsObject() )
    {
      const auto& pd = cobj["port_directions"];
      for ( auto it = pd.MemberBegin(); it != pd.MemberEnd(); ++it )
      {
        if ( it->value.IsString() )
          cell.port_dirs.emplace( it->name.GetString(), it->value.GetString() );
      }
    }

    if ( cobj.HasMember( "parameters" ) && cobj["parameters"].IsObject() )
    {
      const auto& pv = cobj["parameters"];
      for ( auto it = pv.MemberBegin(); it != pv.MemberEnd(); ++it )
      {
        cell.parameters.emplace_back( it->name.GetString(), raw_to_string_( it->value ) );
      }
    }
    if ( cobj.HasMember( "attributes" ) && cobj["attributes"].IsObject() )
    {
      const auto& av = cobj["attributes"];
      for ( auto it = av.MemberBegin(); it != av.MemberEnd(); ++it )
      {
        cell.attributes.emplace_back( it->name.GetString(), raw_to_string_( it->value ) );
      }
    }

    if ( cobj.HasMember( "connections" ) && cobj["connections"].IsObject() )
    {
      const auto& conn = cobj["connections"];
      for ( auto it = conn.MemberBegin(); it != conn.MemberEnd(); ++it )
      {
        cell.connections.emplace( it->name.GetString(), read_bits_( it->value ) );
      }
    }

    ++c_it_;
    out = std::move( cell );
    return instance_return_code::valid;
  }

  instance_return_code read_net_( instance_t& out )
  {
    const auto& mem = *n_it_;
    const std::string name = mem.name.GetString();
    const auto& nobj = mem.value;
    if ( !nobj.IsObject() )
    {
      ++n_it_;
      return instance_return_code::invalid;
    }

    net_name_instance_t net;
    net.name = name;
    if ( nobj.HasMember( "bits" ) )
      net.bits = read_bits_( nobj["bits"] );

    ++n_it_;
    out = std::move( net );
    return instance_return_code::valid;
  }

private:
  bool parse_ok_ = true;
  std::string module_name_;

  rapidjson::Document doc_;
  const rapidjson::Value* module_ = nullptr;

  const rapidjson::Value* ports_ = nullptr;
  const rapidjson::Value* cells_ = nullptr;
  const rapidjson::Value* nets_ = nullptr;
  lorina::diagnostic_engine* diag_;

  rapidjson::Value::ConstMemberIterator p_it_{}, p_end_{};
  rapidjson::Value::ConstMemberIterator c_it_{}, c_end_{};
  rapidjson::Value::ConstMemberIterator n_it_{}, n_end_{};
};

} // namespace json

} // namespace io

} // namespace rinox