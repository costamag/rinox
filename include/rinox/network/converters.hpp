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
    if constexpr ( DesignType == design_type_t::CELL_BASED )
      children.push_back( old_to_new[ni_src][fi_src.output] );
    else if constexpr ( DesignType == design_type_t::ARRAY_BASED )
      children.push_back( old_to_new[ni_src][0] );
    else
      throw std::runtime_error( "[ERROR] invalid design type" );
  } );

  std::vector<unsigned int> ids;
  if constexpr ( DesignType == design_type_t::CELL_BASED )
  {
    auto cell = ntk_src.get_cell( n_src );
    for ( auto const& g : cell.gates )
      ids.push_back( g.id );
  }
  else if constexpr ( DesignType == design_type_t::ARRAY_BASED )
  {
    auto tt = ntk_src.node_function( n_src );
    //auto shrinked = kitty::shrink_to( tt, children.size() );
    ids = { ntk.insert( tt ) };
  }
  else
    throw std::runtime_error( "[ERROR] invalid design type" );

  auto f = ntk.create_node( children, ids );
  std::vector<typename BndNtk::signal> fs;
  for ( auto i = 0u; i < ids.size(); ++i )
    fs.push_back( typename BndNtk::signal{ f.index, i } );
  old_to_new[n_src] = fs;
  ntk_src.set_value( n_src, ntk_src.trav_id() );
}

template<class Ntk, design_type_t DesignType = design_type_t::CELL_BASED, uint32_t MaxNumOutputs = 2u>
void convert_mapped_to_bound( bound_network<DesignType, MaxNumOutputs>& ntk, Ntk const& ntk_src )
{
  using node = typename bound_network<DesignType, MaxNumOutputs>::node;
  using signal = typename bound_network<DesignType, MaxNumOutputs>::signal;
  mockturtle::unordered_node_map<std::vector<signal>, Ntk> old_to_new( ntk_src );

  /* create te constants */
  old_to_new[ntk_src.get_node( ntk_src.get_constant( false ) )] = { ntk.get_constant( false ) };
  old_to_new[ntk_src.get_node( ntk_src.get_constant( true ) )] = { ntk.get_constant( true ) };
  /* create the pis */
  ntk_src.incr_trav_id();
  ntk_src.foreach_pi( [&]( auto f_src ) {
    auto const n_src = ntk_src.get_node( f_src );
    auto const f = ntk.create_pi();
    old_to_new[n_src] = { f };
    ntk_src.set_value( n_src, ntk_src.trav_id() );
  } );

  kitty::dynamic_truth_table buf( 1u );
  kitty::dynamic_truth_table inv( 1u );
  kitty::create_from_binary_string( buf, "10" );
  kitty::create_from_binary_string( inv, "01" );
  /* create the other nodes */

  ntk_src.foreach_po( [&]( auto fo_src ) {
    convert_mapped_to_bound_tfi( ntk, ntk_src, old_to_new, fo_src );
    auto no_src = ntk_src.get_node( fo_src );
    if constexpr ( DesignType == design_type_t::CELL_BASED )
      ntk.create_po( old_to_new[no_src][fo_src.output] );
    else if constexpr ( DesignType == design_type_t::ARRAY_BASED )
      ntk.create_po( old_to_new[no_src][0] );
    else
      throw std::runtime_error( "[ERROR] invalid design type" );
  } );
}

template<class Ntk>
void convert_bound_to_klut(mockturtle::klut_network& klut, Ntk const& ntk)
{
  using node        = typename Ntk::node;
  using klut_signal = mockturtle::klut_network::signal;

  // Signal -> klut signal map (no implicit default-insert!)
  rinox::network::incomplete_signal_map<klut_signal, Ntk> old_to_new(ntk);

  auto fanins_all_mapped = [&](node n) -> bool {
    bool ok = true;
    ntk.foreach_fanin(n, [&](auto fi, auto){
      ok &= old_to_new.has(fi);
    });
    return ok;
  };

  // 1) Seed constants and PIs
  old_to_new[ ntk.get_constant(false) ] = klut.get_constant(false);
  old_to_new[ ntk.get_constant(true)  ] = klut.get_constant(true);

  ntk.foreach_pi([&](auto n_pi){
    old_to_new[ ntk.make_signal(n_pi) ] = klut.create_pi();
  });

  // 2) Worklist: gates whose fanins are all mapped
  std::queue<node> q;
  std::vector<uint8_t> enqueued(ntk.size(), 0);
  std::vector<uint8_t> done(ntk.size(), 0);

  ntk.foreach_gate([&](auto n){
    if (fanins_all_mapped(n)) {
      q.push(n);
      enqueued[n] = 1;
    }
  });

  uint32_t mapped_gates = 0;

  auto try_enqueue_fanouts = [&](node n) {
    // Only enqueue *gates* whose *all* fanins now map to klut signals
    std::unordered_set<node> uniq;
    ntk.foreach_fanout(n, [&](auto g){
      if (!ntk.is_pi(g) && !ntk.is_constant(g)) uniq.insert(g);
    });
    for (auto g : uniq) {
      if (!done[g] && !enqueued[g] && fanins_all_mapped(g)) {
        q.push(g);
        enqueued[g] = 1;
      }
    }
  };

  // 3) Process the queue
  while (!q.empty()) {
    auto n = q.front(); q.pop();
    if (done[n]) continue;

    // Build children in the given fanin order; refuse to proceed if any missing
    std::vector<klut_signal> children;
    ntk.foreach_fanin(n, [&](auto fi, auto){
      if (!old_to_new.has(fi)) {
        throw std::runtime_error("convert_bound_to_klut: fanin not mapped (cycle or scheduling bug)");
      }
      children.push_back(old_to_new[fi]);
    });

    // Create all outputs of this (possibly multi-output) node
    auto const tts = ntk.get_functions(n);
    ntk.foreach_output(n, [&](auto const& f){
      old_to_new[f] = klut.create_node(children, tts[f.output]);
    });

    done[n] = 1;
    ++mapped_gates;

    // See if this unlocks any fanouts
    try_enqueue_fanouts(n);
  }

  // 4) If not all gates were mapped, there’s a real cycle
  //if (mapped_gates < ntk.num_gates()) {
  //  throw std::runtime_error("Cycle detected in network");
  // /}

  // 5) POs
  ntk.foreach_po([&](auto f){
    if (!old_to_new.has(f))
      throw std::runtime_error("convert_bound_to_klut: PO unmapped");
    klut.create_po(old_to_new[f]);
  });
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
template<
    class Ntk,
    network::design_type_t DesignType,
    std::enable_if_t<DesignType == network::design_type_t::CELL_BASED, int> = 0,
    uint32_t MaxNumOutputs = 2u>
bound_network<DesignType, MaxNumOutputs> convert_mapped_to_bound( Ntk const& ntk_src, std::vector<mockturtle::gate> const& gates )
{
  bound_network<DesignType, MaxNumOutputs> ntk( gates );
  convert_mapped_to_bound( ntk, ntk_src );
  return ntk;
}

template<
    class Ntk,
    network::design_type_t DesignType,
    std::enable_if_t<DesignType == network::design_type_t::ARRAY_BASED, int> = 0,
    uint32_t MaxNumOutputs = 2u>
bound_network<DesignType, MaxNumOutputs> convert_mapped_to_bound( Ntk const& ntk_src )
{
  bound_network<DesignType, MaxNumOutputs> ntk;
  convert_mapped_to_bound( ntk, ntk_src );
  return ntk;
}

mockturtle::klut_network convert_bound_to_klut( network::bound_network<network::design_type_t::ARRAY_BASED, 2> & ntk )
{
  mockturtle::klut_network klut;
  convert_bound_to_klut( klut, ntk );
  return klut;
}

} // namespace rinox::network