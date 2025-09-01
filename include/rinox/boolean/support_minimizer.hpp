#pragma once
#include <algorithm>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <kitty/kitty.hpp>
#include <vector>

namespace rinox
{

namespace boolean
{

// Bron–Kerbosch with pivoting for up to 64-node graphs
template<size_t N = 64>
class max_clique_finder
{
public:
  using bitset = std::bitset<N>;

  void set_graph( size_t num_nodes, const std::vector<std::pair<size_t, size_t>>& edges )
  {
    n_ = num_nodes;
    for ( size_t i = 0; i < n_; ++i )
    {
      adj_[i].reset();
    }
    for ( auto [u, v] : edges )
    {
      adj_[u].set( v );
      adj_[v].set( u );
    }
  }

  std::vector<size_t> find()
  {
    best_size_ = 0;
    best_.reset();
    bitset R, P, X;
    for ( size_t i = 0; i < n_; ++i )
      P.set( i );
    expand( R, P, X );
    std::vector<size_t> result;
    for ( size_t i = 0; i < n_; ++i )
      if ( best_.test( i ) )
        result.push_back( i );
    return result;
  }

private:
  void expand( bitset R, bitset P, bitset X )
  {
    if ( P.none() && X.none() )
    {
      if ( R.count() > best_size_ )
      {
        best_ = R;
        best_size_ = R.count();
      }
      return;
    }

    size_t pivot = 0;
    bitset union_PX = P | X;
    for ( size_t i = 0; i < n_; ++i )
      if ( union_PX.test( i ) )
      {
        pivot = i;
        break;
      }

    bitset candidates = P & ~adj_[pivot];
    for ( size_t i = 0; i < n_; ++i )
    {
      if ( !candidates.test( i ) )
        continue;
      R.set( i );
      expand( R & adj_[i], P & adj_[i], X & adj_[i] );
      R.reset( i );
      P.reset( i );
      X.set( i );
    }
  }

  size_t n_ = 0;
  bitset adj_[N];
  bitset best_;
  size_t best_size_ = 0;
};

template<uint32_t MaxNumVars, bool Exact = false>
class support_minimizer
{
public:
  using trut_t = kitty::static_truth_table<MaxNumVars>;
  using func_t = kitty::ternary_truth_table<trut_t>;

  support_minimizer()
  {
    vars_.reserve( MaxNumVars );
    funcs_.reserve( MaxNumVars );
    in_clique_.reserve( MaxNumVars );
    pos_.resize( 256, invalid_ );
  }

  template<typename Vec, typename... Containers>
  void run( func_t& tt, Vec& support, Containers&... extra )
  {
    const size_t sz = support.size();

    if constexpr ( !Exact )
    {
      auto swap_pop_at = []( auto& v, size_t i ) { if (i+1!=v.size()) std::swap(v[i], v.back()); v.pop_back(); };
      for ( size_t i = support.size(); i-- > 0; )
      {
        const auto x = static_cast<uint32_t>( support[i] );
        const auto care0 = kitty::cofactor0( tt._care, x );
        const auto care1 = kitty::cofactor1( tt._care, x );
        const auto care = care0 & care1;
        const auto bits0 = kitty::cofactor0( tt._bits, x );
        const auto bits1 = kitty::cofactor1( tt._bits, x );
        if ( kitty::equal( bits0 & care, bits1 & care ) )
        {
          tt._bits = ( bits0 & care0 ) | ( bits1 & care1 );
          tt._care = care0 | care1;
          swap_pop_at( support, i );
          ( swap_pop_at( extra, i ), ... );
        }
      }
      return;
    }

    vars_.clear();
    funcs_.clear();
    for ( size_t i = 0; i < sz; ++i )
    {
      const auto x = static_cast<uint32_t>( support[i] );
      const auto care0 = kitty::cofactor0( tt._care, x );
      const auto care1 = kitty::cofactor1( tt._care, x );
      const auto care = care0 & care1;
      const auto bits0 = kitty::cofactor0( tt._bits, x );
      const auto bits1 = kitty::cofactor1( tt._bits, x );
      if ( kitty::equal( bits0 & care, bits1 & care ) )
      {
        trut_t bits = ( bits0 & care0 ) | ( bits1 & care1 );
        trut_t cc = care0 | care1;
        vars_.push_back( support[i] );
        func_t f;
        f._bits = bits;
        f._care = cc;
        funcs_.push_back( std::move( f ) );
      }
    }
    const size_t m = vars_.size();
    if ( !m )
      return;

    max_clique_finder<64> mc;
    std::vector<std::pair<size_t, size_t>> edges;

    for ( size_t i = 0; i + 1 < m; ++i )
    {
      for ( size_t j = i + 1; j < m; ++j )
      {
        const auto cij = funcs_[i]._care & funcs_[j]._care;
        if ( kitty::equal( funcs_[i]._bits & cij, funcs_[j]._bits & cij ) )
          edges.emplace_back( i, j );
      }
    }

    mc.set_graph( m, edges );
    const auto& clique = mc.find();

    trut_t acc_bits;
    trut_t acc_care;
    acc_bits = funcs_[clique[0]]._bits;
    acc_care = funcs_[clique[0]]._care;
    for ( size_t idx = 1; idx < clique.size(); ++idx )
    {
      const auto& g = funcs_[clique[idx]];
      trut_t take_g = g._care & ~acc_care;
      acc_bits = ( acc_bits & acc_care ) | ( g._bits & take_g );
      acc_care = acc_care | g._care;
    }
    tt._bits = acc_bits;
    tt._care = acc_care;

    std::fill( pos_.begin(), pos_.end(), invalid_ );
    for ( size_t i = 0; i < sz; ++i )
      pos_[support[i]] = i;
    in_clique_.assign( sz, 0 );
    for ( auto v : clique )
    {
      const auto var = vars_[static_cast<size_t>( v )];
      const auto p = pos_[var];
      if ( p != invalid_ )
        in_clique_[p] = 1;
    }
    auto swap_pop_at = []( auto& v, size_t i ) { if (i+1!=v.size()) std::swap(v[i], v.back()); v.pop_back(); };
    for ( size_t i = support.size(); i-- > 0; )
      if ( in_clique_[i] )
      {
        swap_pop_at( support, i );
        ( swap_pop_at( extra, i ), ... );
      }
  }

private:
  static constexpr size_t invalid_ = static_cast<size_t>( -1 );
  std::vector<func_t> funcs_;
  std::vector<uint8_t> vars_;
  std::vector<char> in_clique_;
  std::vector<size_t> pos_;
};

} // namespace boolean

} // namespace rinox
