/* mockturtle: C++ logic network library
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
  \file mapped_to_bound.hpp
  \brief Transforms a mapped logic network into a bound network

  \author Andrea Costamagna
*/

#pragma once

#include "network.hpp"

namespace rinox::network
{

template<class Ntk, class Map, design_type_t DesignType = design_type_t::CELL_BASED, uint32_t MaxNumOutputs = 2u>
void convert_mapped_to_bound_tfi( bound_network<DesignType, MaxNumOutputs>& ntk, Ntk& ntk_src, Map& old_to_new, typename Ntk::signal f_src )
{
  using BndNtk = bound_network<DesignType, MaxNumOutputs>;

  typename Ntk::node n_src = ntk_src.get_node( f_src );
  if ( ntk_src.value( n_src ) == ntk_src.trav_id() || ntk_src.is_pi( n_src ) )
  {
    return;
  }

  std::vector<typename BndNtk::signal> children;
  ntk_src.foreach_fanin( n_src, [&]( auto fi_src ) {
    auto const ni_src = ntk_src.get_node( fi_src );
    convert_mapped_to_bound_tfi( ntk, ntk_src, old_to_new, fi_src );
    children.push_back( old_to_new[ni_src][fi_src.output] );
  } );
  auto cell = ntk_src.get_cell( n_src );
  std::vector<unsigned int> ids;
  for ( auto const& g : cell.gates )
    ids.push_back( g.id );

  auto f = ntk.create_node( children, ids );
  std::vector<typename BndNtk::signal> fs;
  for ( auto i = 0u; i < ids.size(); ++i )
    fs.push_back( typename BndNtk::signal{ f.index, i } );
  old_to_new[n_src] = fs;
  ntk_src.set_value( n_src, ntk_src.trav_id() );
}

template<class Ntk, design_type_t DesignType = design_type_t::CELL_BASED, uint32_t MaxNumOutputs = 2u>
void convert_mapped_to_bound( bound_network<DesignType, MaxNumOutputs>& ntk, Ntk const& ntk_src, std::vector<mockturtle::gate> const& gates )
{
  using node = typename bound_network<DesignType, MaxNumOutputs>::node;
  using signal = typename bound_network<DesignType, MaxNumOutputs>::signal;
  mockturtle::unordered_node_map<std::vector<signal>, Ntk> old_to_new( ntk_src );
  /* create the pis */
  ntk_src.incr_trav_id();
  ntk_src.foreach_pi( [&]( auto f_src ) {
    auto const n_src = ntk_src.get_node( f_src );
    auto const f = ntk.create_pi();
    old_to_new[n_src] = { f };
    ntk_src.set_value( n_src, ntk_src.trav_id() );
  } );

  kitty::dynamic_truth_table buf(1u);
  kitty::dynamic_truth_table inv(1u);
  kitty::create_from_binary_string( buf, "10" );
  kitty::create_from_binary_string( inv, "01" );
  /* create the other nodes */
  ntk_src.foreach_po( [&]( auto fo_src ) {
    auto no_src = ntk_src.get_node( fo_src );
//    if ( ntk_src.is_constant( no_src ) )
//    {
//      node n;
//      if ( ntk_src.constant_value( no_src ) )
//        
//      else
//    }
//    else
    {
      convert_mapped_to_bound_tfi( ntk, ntk_src, old_to_new, fo_src );
      ntk.create_po( old_to_new[no_src][fo_src.output] );
    }
  } );
}

/*! \brief Convert a mapped network into a bound nework
 *
 * This function is a wrapper function for resynthesizing a k-LUT network (type `NtkSrc`) into a
 * new graph (of type `NtkDest`). The new data structure can be of type AIG, XAG, MIG or XMG.
 * First the function attempts a Disjoint Support Decomposition (DSD), branching the network into subnetworks.
 * As soon as DSD can no longer be done, there are two possibilities depending on the dimensionality of the
 * subnetwork to be resynthesized. On the one hand, if the size of the associated support is lower or equal
 * than 4, the solution can be recovered by exploiting the mapping of the subnetwork to its NPN-class.
 * On the other hand, if the support size is higher than 4, A Shannon decomposition is performed, branching
 * the network in further subnetworks with reduced support.
 * Finally, once the threshold value of 4 is reached, the NPN mapping completes the graph definition.
 *
 * \tparam NtkDest Type of the destination network. Its base type should be either `aig_network`, `xag_network`, `mig_network`, or `xmg_network`.
 * \tparam NtkSrc Type of the source network. Its base type should be `klut_network`.
 * \param ntk_src Input k-lut network
 * \return An equivalent AIG, XAG, MIG or XMG network
 */
template<class Ntk, design_type_t DesignType = design_type_t::CELL_BASED, uint32_t MaxNumOutputs = 2u>
bound_network<DesignType, MaxNumOutputs> convert_mapped_to_bound( Ntk const& ntk_src, std::vector<mockturtle::gate> const& gates )
{
  bound_network<DesignType, MaxNumOutputs> ntk( gates );
  convert_mapped_to_bound( ntk, ntk_src, gates );
  return ntk;
}

} // namespace rinox::network