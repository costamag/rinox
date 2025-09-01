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
  \file permutation.hpp
  \brief Implements utilities to handle permutations

  \author Andrea Costamagna
*/

#pragma once

#include <kitty/kitty.hpp>

namespace rinox
{

namespace boolean
{

/*! \brief Store permutation transformations for up to 16 inputs */
struct permutation_t
{
  uint64_t fmap = 0;
  uint64_t imap = 0;
  uint8_t num_vars = 0;

  permutation_t() = default;
  permutation_t( permutation_t&& ) = default;
  permutation_t( permutation_t const& ) = default;
  permutation_t( std::vector<uint8_t> const& perm )
  {
    num_vars = perm.size();
    for ( uint8_t i{ 0u }; i < perm.size(); ++i )
    {
      set( i, perm[i] );
    }
  }

  permutation_t& operator=( permutation_t const& ) = default;
  permutation_t& operator=( permutation_t&& ) = default;

  constexpr uint8_t forward( uint8_t i ) const noexcept
  {
    return ( fmap >> ( 4 * i ) ) & 0xF;
  }

  constexpr uint8_t inverse( uint8_t i ) const noexcept
  {
    return ( imap >> ( 4 * i ) ) & 0xF;
  }

  constexpr void set( uint8_t i, uint8_t v ) noexcept
  {
    fmap |= static_cast<uint64_t>( v ) << ( 4 * i );
    imap |= static_cast<uint64_t>( i ) << ( 4 * v );
  }

  constexpr bool operator==( const permutation_t& other ) const noexcept
  {
    return ( fmap == other.fmap ) && ( imap == other.imap );
  }
};

template<typename Vec>
Vec inverse_permute( permutation_t const& perm, Vec const& vec )
{
  assert( vec.size() < 16 );
  Vec res( vec.size() );
  for ( auto i = 0u; i < vec.size(); ++i )
  {
    res[i] = vec[perm.inverse( i )];
  }
  return res;
}

template<typename Vec>
Vec forward_permute( permutation_t const& perm, Vec const& vec )
{
  assert( vec.size() < 16 );
  Vec res;
  res.resize( vec.size() );
  for ( auto i = 0u; i < vec.size(); ++i )
  {
    res[i] = vec[perm.forward( i )];
  }
  return res;
}

template<typename Vec>
void inverse_permute_inplace_single( permutation_t const& perm, Vec& vec )
{
  vec = inverse_permute( perm, vec );
}

template<typename... Vecs>
void inverse_permute_inplace( permutation_t const& perm, Vecs&... vecs )
{
  ( inverse_permute_inplace_single( perm, vecs ), ... );
}

template<typename Vec>
void forward_permute_inplace_single( permutation_t const& perm, Vec& vec )
{
  vec = forward_permute( perm, vec );
}

template<typename... Vecs>
void forward_permute_inplace( permutation_t const& perm, Vecs&... vecs )
{
  ( forward_permute_inplace_single( perm, vecs ), ... );
}

template<typename T, typename Fn>
void foreach_permutation( std::vector<T>& vec, Fn&& fn )
{
  auto const n = vec.size();
  auto perms = [&]( int i, auto&& perms_ref ) -> void {
    if ( i >= n )
    {
      fn( vec );
      return;
    }

    for ( int j = i; j < n; ++j )
    {
      std::swap( vec[i], vec[j] );
      perms_ref( i + 1, perms_ref );
      std::swap( vec[i], vec[j] ); // backtrack
    }
  };

  perms( 0, perms );
}

} // namespace boolean

} // namespace rinox