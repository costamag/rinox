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
  \file exact_lib_manager.hpp
  \brief Manager encapsulating operations on network-based databases.

  \author Andrea Costamagna
*/

#pragma once

#include "../evaluation/chains.hpp"
#include "../evaluation/chains/bound_chain.hpp"
#include "../network/utils.hpp"
#include "../traits.hpp"
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/network_utils.hpp>
#include <mockturtle/utils/struct_library.hpp>
#include <mockturtle/utils/tech_library.hpp>

#include <functional>
#include <optional>
#include <type_traits>

namespace rinox
{

namespace databases
{

namespace dispatch
{
/*! \brief Dispatch the resynthesis engine.
 *
 * **Supported network types:**
 * - `aig_network`
 * - `mig_network`
 * - `xag_network`
 */
template<typename Ntk>
struct resynthesis;

/*! \brief Dispatch the element type: signal or literal.
 *
 * **Supported network types:**
 * - `aig_network`
 * - `mig_network`
 * - `xag_network`
 * - `mig_chain`
 * - `xag_chain<true>`
 * - `xag_chain<false>`
 */
template<typename Ntk, typename Enable = void>
struct element;

/*! \brief Dispatch the inversion operation.
 *
 *  Lists and networks use different conventions for the inversion operation.
 *
 * **Supported network types:**
 * - `aig_network`
 * - `mig_network`
 * - `xag_network`
 * - `mig_chain`
 * - `xag_chain<true>`
 * - `xag_chain<false>`
 */
template<typename Ntk>
typename element<Ntk>::type invert( Ntk& ntk, typename element<Ntk>::type const& f )
{
  if constexpr ( mockturtle::is_network_type_v<Ntk> )
    return ntk.create_not( f );
  else if constexpr ( traits::is_boolean_chain<Ntk>::value )
    return ntk.add_not( f );
  else
    static_assert( traits::dependent_false<Ntk>::value,
                   "Unsupported Ntk type in dispatch::invert." );
}

/*! \brief Dispatch the output creation operation.
 *
 *  Lists and networks use different conventions for creating an output.
 *
 * **Supported network types:**
 * - `aig_network`
 * - `mig_network`
 * - `xag_network`
 * - `mig_chain`
 * - `xag_chain<true>`
 * - `xag_chain<false>`
 */
template<typename Ntk>
typename element<Ntk>::type output_inversion( Ntk& ntk,
                                              typename element<Ntk>::type const& f,
                                              bool phase )
{
  auto const res = phase ? invert( ntk, f ) : f;
  return res;
}

template<typename NtkDb, typename NtkDest>
typename element<NtkDest>::type create_node( NtkDest& ntk_dest,
                                             std::vector<typename dispatch::element<NtkDest>::type> const& children,
                                             NtkDb& ntk_db,
                                             mockturtle::node<NtkDb> const& n )
{
  /* XAG */
  if constexpr ( std::is_same<NtkDest, mockturtle::xag_network>::value )
  {
    if ( ntk_db.is_and( n ) )
      return ntk_dest.create_and( children[0], children[1] );
    else
      return ntk_dest.create_xor( children[0], children[1] );
  } /* AIG */
  else if constexpr ( std::is_same<NtkDest, mockturtle::aig_network>::value )
  {
    return ntk_dest.create_and( children[0], children[1] );
  } /* XAG list */
  else if constexpr ( std::is_same<NtkDest, evaluation::chains::xag_chain<true>>::value ||
                      std::is_same<NtkDest, evaluation::chains::xag_chain<false>>::value )
  {
    if ( ntk_db.is_and( n ) )
      return ntk_dest.add_and( children[0], children[1] );
    else
      return ntk_dest.add_xor( children[0], children[1] );
  } /* MIG */
  else if constexpr ( std::is_same<NtkDest, mockturtle::mig_network>::value )
  {
    return ntk_dest.create_maj( children[0], children[1], children[2] );
  } /* MIG list */
  else if constexpr ( std::is_same<NtkDest, evaluation::chains::mig_chain>::value )
  {
    return ntk_dest.add_maj( children[0], children[1], children[2] );
  }
  else
    static_assert( traits::dependent_false<NtkDest>::value,
                   "Unsupported NtkDest type in dispatch::create_node." );
}

template<>
struct resynthesis<mockturtle::aig_network>
{
  using type = mockturtle::xag_npn_resynthesis<mockturtle::aig_network, mockturtle::aig_network, mockturtle::xag_npn_db_kind::aig_complete>;
};

template<>
struct resynthesis<mockturtle::xag_network>
{
  using type = mockturtle::xag_npn_resynthesis<mockturtle::xag_network, mockturtle::xag_network, mockturtle::xag_npn_db_kind::xag_complete>;
};

template<>
struct resynthesis<mockturtle::mig_network>
{
  using type = mockturtle::mig_npn_resynthesis;
};

template<typename Ntk>
struct element<Ntk, std::enable_if_t<mockturtle::is_network_type_v<Ntk>>>
{
  using type = typename Ntk::signal;
};

template<typename Ntk>
struct element<Ntk, std::enable_if_t<traits::is_boolean_chain_v<Ntk>>>
{
  using type = typename Ntk::element_type;
};

}; // namespace dispatch

template<bool UseDCs = true>
struct database_manager_params
{
  /* area of a gate */
  static constexpr float area_gate{ 1.0f };
  /* area of an inverter */
  static constexpr float area_inverter{ 0.0f };
  /* delay of a gate */
  static constexpr float delay_gate{ 1.0f };
  /* delay of an inverter */
  static constexpr float delay_inverter{ 0.0f };
  /* classify in NP instead of NPN */
  static constexpr bool np_classification{ false };
  /* Compute DC classes for matching with  don't cares */
  static constexpr bool compute_dc_classes{ UseDCs };
  /* verbose */
  static constexpr bool verbose{ false };

  static constexpr mockturtle::exact_library_params lib_ps{ area_gate,
                                                            area_inverter,
                                                            delay_gate,
                                                            delay_inverter,
                                                            np_classification,
                                                            compute_dc_classes,
                                                            verbose };
};

/*! \brief Database manager for sub-network reuse.
 *
 * This engine encapsulates operations on network-based databases, i.e.,
 * databases that store sub-networks within a larger network. Given an
 * incompletely specified Boolean function, the `lookup` function identifies
 * matching sub-networks in the database that implement the desired
 * functionality. It returns a unique match type, which contains information
 * about input/output negations and input permutations.
 *
 * Using the returned matching information, the `insert` function reconstructs
 * the matched sub-network in a destination network. The database manager takes
 * care of handling input permutations and inversions, abstracting away the
 * low-level details.
 *
 * \tparam Ntk Database type.
 *
 * \verbatim embed:rst

   Example

   .. code-block:: c++

      static constexpr uint32_t num_vars = 4u;
      using TT = kitty::static_truth_table<num_vars>;
      TT onset = ...;
      database_manager<Ntk> mng;
      kitty::ternary_truth_table<TT> tt( onset );
      auto info = mng.lookup_npn( tt );
      if( info )
      {
        auto sig_out = mng.insert( *info, ntk, f, pis.begin(), pis.end() );
        auto lit_out = mng.insert( *info, list, f, pis_list.begin(), pis_list.end() );
      }

   \endverbatim
 */
template<typename NtkDb, bool UseDCs = false>
class database_manager
{
  static constexpr unsigned num_vars = 4u;
  using library_t = mockturtle::exact_library<NtkDb, num_vars>;
  using supergates_list_t = std::vector<mockturtle::exact_supergate<NtkDb, num_vars>>;
  using node = typename NtkDb::node;
  using signal = typename NtkDb::signal;
  using resyn_t = typename dispatch::resynthesis<NtkDb>::type;

  /*! \brief Type to store the result of Boolean matching.
   */
  struct matches_t
  {
    /*! \brief Iterator over the matching sub-structures.
     *
     *  Each matching sub-structure is uniquely identified by the `root` signal.
     */
    template<typename Fn>
    void foreach_entry( Fn&& fn ) const
    {
      static_assert( std::is_invocable_v<Fn, signal>,
                     "Fn must be callable with signal" );
      for ( auto const& dag : *structures )
      {
        fn( dag.root );
      }
    }

    /*! \brief Encodes the input negations to apply. */
    uint32_t negation;
    /*! \brief Contains the input permutations to apply. */
    std::array<uint8_t, num_vars> permutation;
    /*! \brief Output negation: negate when `phase = true`. */
    bool phase;
    /*! \brief Pointer to a list to identify the sub-networks. */
    supergates_list_t const* structures;
  };

public:
  explicit database_manager()
      : library( *get_cached_library( resyn, ps.lib_ps ) ),
        database( library.get_database() )
  {
  }

  static auto& get_cached_library( const resyn_t& resyn,
                                   const mockturtle::exact_library_params& lib_ps )
  {
    static auto cached_library = std::make_shared<library_t>( resyn, lib_ps );
    return cached_library;
  }

  /*! \brief Boolean matching from an incompletely-specified Boolean function.
   *
   * Perform NPN canonization and Boolean matching with don't cares. The results
   * are wrapped in a `matches_t` object.
   *
   * The input leaves must be provided as a continuous container
   * (e.g., `std::array` or `std::vector`) and accessed via a pair of `begin`
   * and `end` iterators. Permutation and inversion information is provided
   * via a `matches_t` object, typically obtained from a `lookup_*` function.
   *
   * Two callable objects must also be provided:
   * - `invert`: A function that applies inversion to an `Element`.
   * - `get_null`: A function that returns a null or dummy element to fill unused positions.
   *
   * \tparam TT Truth table type to represent onset and careset.
   *
   * \param tt Incompletely-specified Boolean function ( `ternary_truth_table<TT>` )
   *
   * \return A `matches_t` object containing the matching sub-networks.
   */
  template<typename TT>
  [[nodiscard]] std::optional<matches_t> lookup_npn( kitty::ternary_truth_table<TT> const& tt )
  {
    matches_t info;
    /* perform NPN canonization of the on-set and don't care-set. */
    auto [tt_npn, neg, perm] = kitty::exact_npn_canonization( tt._bits );
    auto const dc_npn = apply_npn_transformation( ~tt._care, neg & ~( 1 << num_vars ), perm );

    /* extract the sub-network stored as `supergates`. */
    info.structures = library.get_supergates( tt_npn, dc_npn, neg, perm );
    if ( info.structures == nullptr )
    {
      return std::nullopt;
    }

    /* save the input negations and permutations to apply. */
    info.negation = 0;
    for ( auto j = 0u; j < num_vars; ++j )
    {
      info.permutation[perm[j]] = j;
      info.negation |= ( ( neg >> perm[j] ) & 1 ) << j;
    }

    /* save the output negation to apply */
    info.phase = ( ( neg >> num_vars ) & 1 );
    return info;
  }

  /*! \brief Permute and invert leaves to match a target Boolean function.
   *
   * This function applies a permutation and optional inversion to a list of
   * leaves (typically primary inputs), to match the canonical form of a
   * Boolean function.
   *
   * The input leaves must be provided as a continuous container
   * (e.g., `std::array` or `std::vector`) and accessed via a pair of `begin`
   * and `end` iterators. Permutation and inversion information is provided via
   * a `matches_t` object, typically obtained from a `lookup_*` function.
   *
   * Two callable objects must also be provided:
   * - `invert`: A function that applies inversion to an `Element`.
   * - `get_null`: A function that returns a null or dummy element to fill
   * unused positions.
   *
   * \tparam Element Type of the elements to permute/invert (e.g., signals or literals).
   * \tparam Iterator Iterator type for the input container.
   * \tparam FnInv Callable type for inversion.
   * \tparam FnNull Callable type for null element generation.
   *
   * \param info Matching information (e.g., permutation and inversion).
   * \param begin Iterator to the beginning of the input leaves.
   * \param end Iterator to the end of the input leaves.
   * \param invert Function to apply inversion.
   * \param get_null Function to return a null element.
   *
   * \return A `std::array` of size `num_vars` with permuted and conditionally inverted elements.
   */
  template<typename Element, typename Iterator, typename FnInv, typename FnNull>
  std::array<Element, num_vars> match_leaves( matches_t const& info,
                                              Iterator begin,
                                              Iterator end,
                                              FnInv&& invert,
                                              FnNull&& get_null )
  {
    size_t const num_inputs = end - begin;
    assert( num_inputs <= num_vars );
    std::array<Element, num_vars> leaves;
    auto i = 0u;

    /* permutation */
    for ( auto it = begin; it != end; ++it )
    {
      leaves[info.permutation[i++]] = *it;
    }

    while ( i < num_vars )
    {
      leaves[info.permutation[i++]] = get_null();
    }

    /* inversion */
    for ( i = 0u; i < num_vars; ++i )
    {
      if ( ( info.negation >> i ) & 1 )
      {
        leaves[i] = invert( leaves[i] );
      }
    }

    return leaves;
  }

  /*! \brief Insert a database sub-network in a destination network.
   *
   * This function applies matching transformations to a sub-network stored in
   * the database to insert it into a destination network, given the nodes of
   * this network that serve as the inputs of the stored sub-network, upon
   * permutation and inversion.
   *
   * \tparam NtkDest Network type where to store the sub-network.
   * \tparam Iterator Iterator type for the input container.
   *
   * \param matches_t Matching information (e.g., permutation and inversion).
   * \param ntk_dest Network where to store the sub-network.
   * \param root NtkDb signal to extract the desired sub-network from the database.
   * \param begin Iterator to the first leaf.
   * \param end Iterator to the last leaf.
   *
   * \return A `signal<NtkDest>` when NtkDest is a network, `uint32_t` when NtkDest is a list.
   */
  template<typename NtkDest, typename Iterator>
  auto insert( matches_t const& info,
               NtkDest& ntk_dest,
               signal root,
               Iterator begin,
               Iterator end )
  {
    using element = typename dispatch::element<NtkDest>::type;

    auto cond_invert = [&ntk_dest]( auto const& f, bool complemented ) {
      return complemented ? dispatch::invert( ntk_dest, f ) : f;
    };

    /* permute and invert the leaves */
    auto const leaves = match_leaves<element>(
        info,
        begin,
        end,
        [&ntk_dest]( element const& f ) {
          return dispatch::invert<NtkDest>( ntk_dest, f );
        },
        [&ntk_dest]() {
          return ntk_dest.get_constant( false );
        } );

    /* initialize the hash map connecting database nodes to signals in the destination network. */
    auto map = create_map( ntk_dest, leaves );
    database.incr_trav_id();

    for ( auto i = 0u; i < num_vars; ++i )
    {
      node const n = database.pi_at( i );
      database.set_visited( n, database.trav_id() );
    }

    /* insert the database entry in a list or another network */
    std::function<element( signal const& )> synthesize = [&]( signal const& f ) {
      node n = database.get_node( f );
      if ( database.is_constant( n ) || ( database.visited( n ) == database.trav_id() ) )
      {
        return cond_invert( map[n], database.is_complemented( f ) );
      }
      database.set_visited( n, database.trav_id() );

      std::vector<element> children;
      database.foreach_fanin( n, [&]( auto fi ) {
        element s = synthesize( fi );
        children.push_back( s );
      } );

      map[n] = dispatch::create_node<NtkDb, NtkDest>( ntk_dest, children, database, n );

      return cond_invert( map[n], database.is_complemented( f ) );
    };
    element s = synthesize( root );

    return dispatch::output_inversion( ntk_dest, s, info.phase );
  }

  /*! \brief Evaluate the number of gates in the transitive fanin
   */
  size_t get_cost( signal const& f )
  {
    return network::count_nodes( database, f );
  }

private:
  template<typename NtkDest>
  std::unordered_map<node, typename dispatch::element<NtkDest>::type>
  create_map( NtkDest& ntk_dest,
              std::array<typename dispatch::element<NtkDest>::type, num_vars> const& leaves )
  {
    std::unordered_map<node, typename dispatch::element<NtkDest>::type> map;
    for ( auto i = 0u; i < num_vars; ++i )
    {
      node n = database.pi_at( i );
      map[n] = leaves[i];
    }
    map[database.get_node( database.get_constant( false ) )] = ntk_dest.get_constant( false );
    return map;
  }

private:
  /*! \brief Parameters: for now no need to specify them */
  database_manager_params<UseDCs> ps;
  /*! \brief Resynthesis engine */
  resyn_t resyn;
  /*! \brief Exact library defined for the NtkDb type */
  library_t library;
  /*! \brief Database represented as a network */
  NtkDb& database;
};

} /* namespace databases */

} /* namespace rinox */