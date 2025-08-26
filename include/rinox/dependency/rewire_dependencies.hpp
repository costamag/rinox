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
  \file rewire_dependencies.hpp
  \brief Compute dependencies allowing to rewire the pivot node.

  \author Andrea Costamagna
*/

#pragma once

#include "dependency_cut.hpp"

namespace rinox
{

namespace dependency
{

struct default_rewire_params
{
  static constexpr uint32_t max_cuts_size = 6u;
};

template<class Ntk, typename StaticParams = default_rewire_params>
class rewire_dependencies
{
public:
  static constexpr uint32_t max_cuts_size = StaticParams::max_cuts_size;
  using signal_t = typename Ntk::signal;
  using node_index_t = typename Ntk::node;
  using truth_table_t = kitty::static_truth_table<max_cuts_size>;
  using functionality_t = kitty::ternary_truth_table<truth_table_t>;
  using cut_t = dependency_cut_t<Ntk, max_cuts_size>;

public:
  rewire_dependencies( Ntk& ntk )
      : ntk_( ntk )
  {
  }

  template<typename WinMng, typename WinSim>
  void run( WinMng& window, WinSim& simulator )
  {
    using signature_t = typename WinSim::signature_t;
    window.mark_contained();
    cuts.clear();
    auto const n = window.get_pivot();

    //std::cout << "==========================================================" << std::endl;
    //std::cout << "=                       " << n <<  "                     =" << std::endl;
    //std::cout << "==========================================================" << std::endl;

    if ( ntk_.fanin_size( n ) > max_cuts_size )
      return;

    std::vector<signature_t const*> sim_ptrs;
    std::vector<signal_t> leaves_curr;

    ntk_.foreach_fanin( n, [&]( auto fi, auto ii ) {
      leaves_curr.push_back( fi );
    } );

    std::vector<signature_t> tts_curr;
    ntk_.foreach_output( n, [&]( auto const& f ) {
      tts_curr.push_back( simulator.get( f ) );
    } );

    // Current
    //ntk_.foreach_fanin( n, [&]( auto fii, auto iii ) {
    //  kitty::print_binary( simulator.get( fii ) ); std::cout << " +" << iii << std::endl;
    //} );
    //std::cout << n << " " << ntk_.get_expression( typename Ntk::signal{ n, 0 } ) << std::endl;
    //for ( auto tt : tts_curr )
    //{
    //  kitty::print_binary( tt ); std::cout << " curr" << std::endl;
    //}

    signature_t obs_care = simulator.get_careset();
    ntk_.foreach_fanin( n, [&]( auto fi, auto ii ) {
      signature_t flipped = ~simulator.get( fi );

      sim_ptrs.clear();
      ntk_.foreach_fanin( n, [&]( auto fii, auto iii ) {
        if ( ii == iii )
        {
          sim_ptrs.push_back( &flipped );
      //    kitty::print_binary( *sim_ptrs.back() ); std::cout << " -" << iii << std::endl;
        }
        else
        {
          sim_ptrs.push_back( &simulator.get( fii ) ); 
      //    kitty::print_binary( *sim_ptrs.back() ); std::cout << " +" << iii << std::endl;
        }
      } );

      auto tts_flip = ntk_.compute( n, sim_ptrs );
      //for ( auto tt : tts_flip )
      //{
      //  kitty::print_binary( tt ); std::cout << std::endl;
      //}

      signature_t dc = ~obs_care;
      for ( auto i = 0u; i < tts_curr.size(); ++i )
      {
        dc |= ~( tts_flip[i] ^ tts_curr[i] );
      }
      auto const care = ~dc;
      std::vector<signal_t> leaves = leaves_curr;
      window.foreach_divisor( [&]( auto const& f, auto i ) {
        if ( f != fi )
        {
          auto const& sim_curr = simulator.get( fi );
          auto const& sim_cand = simulator.get( f );
          if ( kitty::equal( sim_curr & care, sim_cand & care ) )
          {
            leaves[ii] = f;
            auto func = extract_function<signature_t, max_cuts_size>( sim_ptrs, sim_curr, care );
            cuts.emplace_back( dependency_t::REWIRE_DEP, n, leaves, func );
          }
        }
      } );

      sim_ptrs[ii] = &simulator.get( fi );
    } );
  }

  template<typename Fn>
  void foreach_cut( Fn&& fn )
  {
    for ( auto i = 0u; i < cuts.size(); ++i )
    {
      fn( cuts[i], i );
    }
  }

private:
  Ntk& ntk_;
  std::vector<cut_t> cuts;
};

} // namespace dependency

} // namespace rinox