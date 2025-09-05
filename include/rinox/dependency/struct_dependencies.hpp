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
  \file struct_dependencies.hpp
  \brief Compute dependencies allowing to window the pivot node.

  \author Andrea Costamagna
*/

#pragma once

#include "../math/math.hpp"
#include "dependency_cut.hpp"

namespace rinox
{

namespace dependency
{

struct default_struct_params
{
  static constexpr uint32_t num_vars_sign = 6u;
  static constexpr uint32_t max_cuts_size = 6u;
  static constexpr uint32_t max_reconv = 3u;
};

template<class Ntk, typename StaticParams = default_struct_params>
class struct_dependencies
{

public:
  static constexpr uint32_t num_vars_sign = StaticParams::num_vars_sign;
  static constexpr uint32_t max_cuts_size = StaticParams::max_cuts_size;
  static constexpr uint32_t max_reconv = StaticParams::max_reconv;

  using signal_t = typename Ntk::signal;
  using signals_t = rinox::math::static_vector<signal_t, max_cuts_size>;
  using node_index_t = typename Ntk::node;
  using signature_t = kitty::static_truth_table<num_vars_sign>;

public:
  struct_dependencies( Ntk& ntk )
      : ntk_( ntk )
  {
  }

  template<typename WinMng, typename WinSim>
  void run( WinMng& window, WinSim& simulator )
  {
    assert( ( std::is_same<signature_t, typename WinSim::signature_t>::value && "signatures have different type" ) );
    cuts_.clear();

    // identify candidates through branch and bound
    struct_enumeration<WinMng, WinSim>( window, simulator );
  }

  template<typename Fn>
  void foreach_cut( Fn&& fn )
  {
    for ( auto i = 0u; i < cuts_.size(); ++i )
    {
      fn( cuts_[i], i );
    }
  }

#pragma region Candidates Identification

  bool contains( signals_t const& cut1, signals_t const& cut2 )
  {
    if ( cut2.size() > cut1.size() )
      return false;
    bool res = true;
    auto i = 0;
    auto j = 0;
    while ( res && ( i < cut1.size() ) && ( j < cut2.size() ) )
    {
      if ( cut1[i] < cut2[j] )
        i++;
      else if ( cut1[i] == cut2[j] )
      {
        i++;
        j++;
      }
      else if ( cut1[i] > cut2[j] )
        res = false;
    }
    if ( res && ( j == cut2.size() ) )
      return true;
    return false;
  }

  bool contains( signals_t const& leaves, std::vector<signals_t> const& leaves_vec )
  {
    for ( auto const& other : leaves_vec )
    {
      if ( ( leaves != other ) && contains( leaves, other ) )
        return true;
    }
    return false;
  }

  // check if cut2 ⊆ cut1 (both sorted)
  bool is_subset( signals_t const& cut1, signals_t const& cut2 )
  {
    if ( cut2.size() > cut1.size() )
      return false;

    size_t i = 0, j = 0;
    while ( i < cut1.size() && j < cut2.size() )
    {
      if ( cut1[i] < cut2[j] )
        ++i;
      else if ( cut1[i] == cut2[j] )
      {
        ++i;
        ++j;
      }
      else
        return false;
    }
    return j == cut2.size();
  }

  template<typename WinMng, typename WinSim>
  void struct_enumeration( WinMng const& window, WinSim& simulator )
  {
    using cut_sigs_t = rinox::math::static_vector<signal_t, max_cuts_size>;
    using front_sigs_t = rinox::math::static_vector<signal_t, max_cuts_size + 10u>;

    struct CutHash
    {
      size_t operator()( front_sigs_t const& cut ) const noexcept
      {
        size_t h = 1469598103934665603ull; // FNV-1a offset
        for ( auto const& s : cut )
        {
          size_t v = ( static_cast<size_t>( s.index ) << 4 ) ^ s.output;
          h ^= v;
          h *= 1099511628211ull;
        }
        return h;
      }
    };

    struct CutEq
    {
      bool operator()( front_sigs_t const& a, front_sigs_t const& b ) const noexcept
      {
        if ( a.size() != b.size() )
          return false;
        for ( size_t i = 0; i < a.size(); ++i )
          if ( a[i] != b[i] )
            return false;
        return true;
      }
    };

    cuts_.clear();
    auto const care = simulator.get_careset();
    auto const pivot = window.get_pivot();

    // === initial front = pivot fanins ===
    front_sigs_t init_front;
    ntk_.foreach_fanin( pivot, [&]( auto const& fi, auto ) {
      init_front.push_back( fi );
    } );

    std::stable_sort( init_front.begin(), init_front.end() );

    std::vector<front_sigs_t> worklist;
    std::unordered_set<front_sigs_t, CutHash, CutEq> seen;

    worklist.push_back( init_front );
    seen.insert( init_front );

    while ( !worklist.empty() )
    {
      front_sigs_t front = worklist.back();
      worklist.pop_back();

      // === 1) Save this cut if within legal size ===
      if ( front.size() <= max_cuts_size )
      {
        cut_sigs_t cut;
        for ( auto const& s : front )
          cut.push_back( s );

        dependency_cut_t<Ntk, max_cuts_size> dep_cut( dependency_t::STRUCT_DEP, pivot, cut );

        // collect simulator signatures
        std::vector<signature_t const*> in_ptrs;
        for ( auto const& l : cut )
        {
          in_ptrs.push_back( &simulator.get( l ) );
        }

        // collect *all signals* of pivot node (multi-output)
        ntk_.foreach_output( pivot, [&]( auto const& f ) {
          auto func = extract_function<signature_t, max_cuts_size>(
              in_ptrs, simulator.get( f ), care );
          dep_cut.add_func( func );
        } );

        cuts_.push_back( dep_cut );
      }

      // === 2) If frontier too large, prune ===
      if ( front.size() > max_cuts_size + max_reconv )
        continue;

      // === 3) Expand each leaf ===
      for ( uint32_t i = 0; i < front.size(); ++i )
      {
        auto leaf_node = ntk_.get_node( front[i] );
        if ( ntk_.is_pi( leaf_node ) )
          continue;

        front_sigs_t new_front;
        // copy all except i
        for ( uint32_t j = 0; j < front.size(); ++j )
          if ( j != i )
            new_front.push_back( front[j] );

        bool abort = false;
        ntk_.foreach_fanin( leaf_node, [&]( auto const& fi, auto ) {
          if ( abort )
            return;

          auto ni = ntk_.get_node( fi );
          if ( !window.is_contained( ni ) )
          {
            abort = true;
            return;
          }

          // add *all outputs* of ni
          ntk_.foreach_output( ni, [&]( auto const& fo ) {
            auto it = std::lower_bound( new_front.begin(), new_front.end(), fo );
            if ( it == new_front.end() || *it != fo )
              new_front.insert( it, fo );
          } );
        } );

        if ( abort )
          continue;

        // === 4) Deduplication using hash set ===
        if ( seen.find( new_front ) == seen.end() )
        {
          seen.insert( new_front );
          worklist.push_back( new_front );
        }
      }
    }
  }

#pragma endregion

private:
  Ntk& ntk_;
  std::vector<dependency_cut_t<Ntk, max_cuts_size>> cuts_;
};

} // namespace dependency

} // namespace rinox