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
  \file homo_xpx_analyzer.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../boolean/simd.hpp"
#include "../math/math.hpp"
#include <kitty/kitty.hpp>

namespace rinox
{

namespace dependency
{

enum dependency_t
{
  REWIRE_DEP, // maintain the gate, non-local rewiring
  STRUCT_DEP, // structural dependency
  WINDOW_DEP, // non-structural dependency. Does not require verification
  SIMULA_DEP  // non-structural dependency. Requires verification
};

template<class Ntk, uint32_t MaxCutSize>
struct dependency_cut_t
{

public:
  using signal_t = typename Ntk::signal;
  using signals_t = rinox::math::static_vector<signal_t, MaxCutSize>;
  using node_index_t = typename Ntk::node;
  using truth_table_t = kitty::static_truth_table<MaxCutSize>;
  using functionality_t = kitty::ternary_truth_table<truth_table_t>;

  dependency_cut_t() = default;

  dependency_cut_t( dependency_t const& type, node_index_t const& root, signals_t const& leaves, functionality_t const& func )
      : type( type ), root( root ), leaves( leaves ), func( { func } )
  {}

  dependency_cut_t( dependency_t const& type, node_index_t const& root, signals_t const& leaves )
      : type( type ), root( root ), leaves( leaves )
  {}

  void add_leaf( signal_t const& f )
  {
    leaves.push_back( f );
  }

  void add_func( functionality_t const& tt )
  {
    func.push_back( tt );
  }

  typename signals_t::const_iterator cbegin()
  {
    return leaves.cbegin();
  }

  typename signals_t::const_iterator cend()
  {
    return leaves.cend();
  }

  typename signals_t::iterator begin()
  {
    return leaves.begin();
  }

  typename signals_t::iterator end()
  {
    return leaves.end();
  }

  size_t size() const
  {
    return leaves.size();
  }

  dependency_t type;
  node_index_t root;
  std::vector<functionality_t> func;
  signals_t leaves;
};

template<uint32_t NumVars>
const std::array<kitty::static_truth_table<NumVars>, NumVars>& get_projection_functions()
{
  static const std::array<kitty::static_truth_table<NumVars>, NumVars> arr = [] {
    std::array<kitty::static_truth_table<NumVars>, NumVars> v;
    // Expensive computation here
    for ( auto i = 0u; i < NumVars; ++i )
      kitty::create_nth_var( v[i], i );
    return v;
  }();
  return arr;
}

template<typename Signature, uint32_t NumVars>
kitty::ternary_truth_table<kitty::static_truth_table<NumVars>> extract_function( std::vector<Signature const*> const& sim_ptrs, Signature const& func, Signature const& care )
{
  using truth_table_t = kitty::static_truth_table<NumVars>;

  uint32_t num_minterms = 1u << sim_ptrs.size();
  truth_table_t onset;
  truth_table_t careset;
  auto const& proj_fns = get_projection_functions<NumVars>();
  for ( auto m = 0u; m < num_minterms; ++m )
  {
    Signature tmp_sig;
    truth_table_t tmp_fun;
    tmp_fun = ~tmp_fun;
    tmp_sig = ~tmp_sig;
    for ( auto v = 0u; v < sim_ptrs.size(); ++v )
    {
      bool value = ( ( m >> v ) & 0x1 ) > 0;
      tmp_sig &= value ? *sim_ptrs[v] : ~( *sim_ptrs[v] );
      tmp_fun &= value ? proj_fns[v] : ~( proj_fns[v] );
    }
    if ( kitty::count_ones( care & tmp_sig ) > 0 )
    {
      careset |= tmp_fun;
      if ( kitty::count_ones( care & func & tmp_sig ) > 0 )
        onset |= tmp_fun;
    }
  }
  kitty::ternary_truth_table<truth_table_t> tt( onset, careset );

  return tt;
}

template<uint32_t NumVars>
class function_enumerator
{
public:
  using truth_table_t = kitty::static_truth_table<NumVars>;
  using functionality_t = kitty::ternary_truth_table<truth_table_t>;

  function_enumerator()
  {
    for ( auto i = 0u; i < NumVars; ++i )
      kitty::create_nth_var( proj_funcs_[i], i );
  }

  template<typename Fn>
  void foreach_dont_care_assignment( functionality_t const& func, uint32_t num_vars, Fn&& fn )
  {
    uint32_t const num_minterms = 1u << num_vars;

    enumerate_dont_cares( func, num_vars );

    truth_table_t const care = func._care;
    truth_table_t const bits = boolean::binary_and( care, func._bits );
    truth_table_t tmp;

    uint64_t num_funcs = 1u << dont_cares_[0].size();
    for ( uint64_t f = 0u; f < num_funcs; ++f )
    {
      truth_table_t tt = bits;
      for ( auto i = 0u; i < dont_cares_[0].size(); ++i )
      {
        if ( ( ( f >> i ) & 0x1 ) > 0 )
        {
          auto const m = dont_cares_[0][i];
          boolean::set_ones( tmp );
          for ( auto v = 0u; v < num_vars; ++v )
          {
            if ( ( ( m >> v ) & 0x1 ) > 0 )
              tmp = boolean::binary_and( tmp, proj_funcs_[v] );
            else
              tmp = boolean::binary_and( tmp, boolean::unary_not( proj_funcs_[v] ) );
          }
          tt = boolean::binary_or( tt, tmp );
        }
      }
      fn( tt );
    }
  }

  template<typename Ntk, typename Fn>
  void foreach_dont_care_assignment( dependency_cut_t<Ntk, NumVars> const& cut, Fn&& fn )
  {
    foreach_dont_care_assignment( cut.func[0], cut.inputs.size(), [&]( auto const& tt ) {
      fn( tt );
    } );
  }

private:
  void enumerate_dont_cares( functionality_t const& tt, uint32_t num_vars )
  {
    uint64_t const num_minterms = 1u << num_vars;
    dont_cares_.clear();
    truth_table_t tmp;
    {

      dont_cares_.emplace_back();
      for ( uint64_t m = 0; m < num_minterms; ++m )
      {
        boolean::set_ones( tmp );
        for ( auto i = 0u; i < num_vars; ++i )
        {
          if ( ( ( m >> i ) & 0x1 ) > 0 )
            tmp = boolean::binary_and( tmp, proj_funcs_[i] );
          else
            tmp = boolean::binary_and( tmp, boolean::unary_not( proj_funcs_[i] ) );
        }
        if ( kitty::count_ones( boolean::binary_and( tmp, tt._care ) ) <= 0 )
          dont_cares_.back().push_back( m );
      }
    }
  }

private:
  std::array<truth_table_t, NumVars> proj_funcs_;
  std::vector<std::vector<uint32_t>> dont_cares_;
};

} // namespace dependency

} // namespace rinox