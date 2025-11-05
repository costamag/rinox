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
  \file partial_specs.hpp
  \brief Synthesis engine from partial specifications

  \author Andrea Costamagna
*/

#ifndef RINOX_SYNTHESIS_FROM_PARTIAL_SPECS_HPP_
#define RINOX_SYNTHESIS_FROM_PARTIAL_SPECS_HPP_

namespace rinox
{

namespace synthesis
{

namespace partial
{

template<typename Signal, uint32_t kCapacity>
class BinaryDecisionTree
{
private:
  struct TreeNode{
    Signal v;
    uint32_t t;
    uint32_t e;
  };

public:
  BinaryDecisionTree()
  {
    TreeNode const0 { 0u, 0u, 0u };
    TreeNode const1 { 1u, 1u, 1u };
    data_.push_back(const0);
    data_.push_back(const1);
  };

  void Reset()
  {
    data_.resize( 2u );
  }

  uint32_t AddNode( Signal const& v, uint32_t t, uint32_t e )
  {
    uint32_t index = data_.size();
    TreeNode n { v, t, e };
    data_.push_back( n );
    return index;
  }

  template<typename Fn>
  void ForeachNode(Fn && fn)
  {
    for ( auto i=2u; i < data_.size(); ++i )
    {
      fn( data_[i], i );
    }
  }

  math::static_vector<TreeNode, kCapacity+2u> data_;
};

template<typename List, typename CFunction, uint32_t kMaxNumInputs>
class Learner
{
public:
  using Signal = typename List::Signal;
  using IFuncion = kitty::ternary_truth_table<CFunction>;
  using Simulations = math::static_vector<CFunction, (1u << kMaxNumInputs)+2u>;

  Learner( List& list )
      : list_( list ) {
    sims_.emplace_back();
    sims_.push_back(~sims_[0]);
  }

  uint32_t ComputeInformation( CFunction const& bits, CFunction const& care ) const
  {
    uint32_t const onset_size = kitty::count_ones( care & bits );
    uint32_t const careset_size = kitty::count_ones( care );
    return onset_size * ( careset_size - onset_size );
  }

  uint32_t Synthesize(std::vector<Signal> const& support, CFunction const& bits, CFunction const& care){
    kitty::ternary_truth_table<CFunction> fun(bits, care);
    std::optional<Signal> best;

    int best_info = ComputeInformation( bits, care );
    for (auto const s : support)
    {
      auto const& ttX = list_.GetSimBits( s );
      auto const care_else = care & ~ttX;
      auto const care_then = care & ttX;
      uint32_t const info = ComputeInformation( bits, care_else ) +
                            ComputeInformation( bits, care_then );

      if (info<=best_info)
      {
        best_info = info;
        best = s;
      }
    }
    auto const& ttX = list_.GetSimBits( *best );
    auto const& care_then = care &  ttX;
    auto const& care_else = care & ~ttX;
    uint32_t t;
    uint32_t e;
    if ( kitty::count_ones( care_then & bits ) == 0 )
      t = 0u;
    else if ( kitty::equal( care_then & bits, care_then ) )
      t = 1u;
    else
      t = Synthesize( support, bits, care_then );

    if ( kitty::count_ones( care_else & bits ) == 0 )
      e = 0u;
    else if ( kitty::equal( care_else & bits, care_else ) )
      e = 1u;
    else
      e = Synthesize( support, bits, care_else );
    
    return bdd_.AddNode(*best, t, e); 
  }

CFunction Run(CFunction const& bits, CFunction const& care, std::vector<Signal> const& support)
{
  bdd_.Reset();
  sims_.resize( 2u );
  Synthesize(support, bits, care);

  CFunction res;
  bdd_.ForeachNode( [&](auto const& n, auto const& index){
    auto const& sim_var = list_.GetSimBits( n.v );
    auto const& sim_then = sims_[n.t];
    auto const& sim_else = sims_[n.e];
    sims_.push_back( (sim_var & sim_then) | (~sim_var & sim_else) );
  } );

  return sims_[sims_.size()-1];
}

private:
  List& list_;
  BinaryDecisionTree<Signal, 1u << kMaxNumInputs> bdd_;
  Simulations sims_;

};

}

}

}

#endif