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
  \file resynthesize.hpp
  \brief Inplace rewrite

  \author Andrea Costamagna
*/

#pragma once

#include "../../databases/mapped_database.hpp"
#include "../../dependency/dependency_cut.hpp"
#include "../../dependency/rewire_dependencies.hpp"
#include "../../dependency/struct_dependencies.hpp"
#include "../../dependency/window_dependencies.hpp"
#include "../../synthesis/lut_decomposer.hpp"
#include "../../windowing/window_manager.hpp"
#include "../../windowing/window_simulator.hpp"
#include "../profilers/area_profiler.hpp"
#include "../profilers/delay_profiler.hpp"
#include "../profilers/power_profiler.hpp"
#include "../profilers/profilers_utils.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <optional>

#ifndef RINOX_NUM_VARS_SIGN
#define RINOX_NUM_VARS_SIGN 6
#endif

#ifndef RINOX_MAX_CUTS_SIZE
#define RINOX_MAX_CUTS_SIZE 6
#endif

#ifndef RINOX_MAX_CUBE_SPFD
#define RINOX_MAX_CUBE_SPFD 12
#endif

#ifndef RINOX_MAX_NUM_LEAVES
#define RINOX_MAX_NUM_LEAVES 12
#endif

namespace rinox
{

namespace opto
{

namespace algorithms
{

/*! \brief Statistics for rewrite.
 *
 * The data structure `rewrite_stats` provides data collected by running
 * `rewrite`.
 */
struct resynthesis_stats
{
  windowing::window_manager_stats window_st;
  /*! \brief Total runtime. */
  mockturtle::stopwatch<>::duration time_total{ 0 };

  /*! \brief Expected gain. */
  uint32_t estimated_gain{ 0 };

  /*! \brief Candidates */
  uint32_t candidates{ 0 };

  uint32_t num_struct{ 0 };
  uint32_t num_window{ 0 };
  uint32_t num_simula{ 0 };
  uint32_t num_rewire{ 0 };
  uint32_t num_optims{ 0 };

  void report() const
  {
    std::cout << fmt::format( "[i] total time       = {:>5.2f} secs\n", mockturtle::to_seconds( time_total ) );
    std::cout << fmt::format( "    num struct       = {:5d}\n", num_struct );
    std::cout << fmt::format( "    num window       = {:5d}\n", num_window );
    std::cout << fmt::format( "    num simula       = {:5d}\n", num_simula );
    std::cout << fmt::format( "    num rewire       = {:5d}\n", num_rewire );
  }
};

template<uint32_t NumLeaves = RINOX_MAX_NUM_LEAVES>
struct default_resynthesis_params
{
  static constexpr uint32_t max_num_leaves = NumLeaves;

  profilers::profiler_params profiler_ps;

  /*! \brief If true, candidates are only accepted if they do not increase logic depth. */
  struct window_manager_params : windowing::default_window_manager_params
  {
    static constexpr uint32_t max_num_leaves = NumLeaves;
    bool preserve_depth = false;
    int32_t odc_levels = 0u;
    uint32_t skip_fanout_limit_for_divisors = 100u;
    uint32_t max_num_divisors{ 128 };
  };

  window_manager_params window_manager_ps;
  uint32_t stop_after_nth{ std::numeric_limits<uint32_t>::max() };

  static constexpr bool do_strashing = true;
  /*! \brief Use satisfiability don't cares for optimization. */
  bool use_dont_cares = false;

  /*! \brief If true try fanin rewiring */
  bool try_rewire = false;

  /*! \brief If true try cut-rewriting with structural cuts */
  bool try_struct = false;

  /*! \brief If true try window-based rewriting with non-structural cuts */
  bool try_window = false;

  /*! \brief If true try simulation-guided rewriting with non-structural cuts */
  bool try_simula = false;

  /*! \brief Activates lazy man's synthesis when set to true */
  bool dynamic_database = false;

  /*! \brief Cube size for the signatures in simulation-guided resubstitution */
  static constexpr uint32_t num_vars_sign = RINOX_NUM_VARS_SIGN;
  /*! \brief Maximum number of leaves in the dependency cuts */
  static constexpr uint32_t max_cuts_size = RINOX_MAX_CUTS_SIZE;
  /*! \brief Maximum cube size exactly represented with SPFDs */
  static constexpr uint32_t max_cube_spfd = RINOX_MAX_CUBE_SPFD;

  /*! \brief Maximum fanout size for a node to be optimized*/
  uint32_t fanout_limit = 12u;
};

namespace detail
{

template<class Ntk, typename Database, typename Profiler, typename Params = default_resynthesis_params<RINOX_MAX_NUM_LEAVES>>
class resynthesize_impl
{
  using node = typename Ntk::node;
  using node_t = typename Ntk::node;
  using signal = typename Ntk::signal;
  using signal_t = typename Ntk::signal;
  using node_index_t = typename Ntk::node;
  using cut_t = dependency::dependency_cut_t<Ntk, Params::max_cuts_size>;
  using chain_t = evaluation::chains::bound_chain<Ntk::design_t>;
  using func_t = kitty::static_truth_table<Params::max_cuts_size>;
  using data_t = kitty::static_truth_table<Database::max_num_vars>;
  using decomposer_t = synthesis::lut_decomposer<Params::max_cuts_size, Database::max_num_vars>;
  using window_manager_t = windowing::window_manager<Ntk, typename Params::window_manager_params>;
  // using validator_t = circuit_validator<Ntk, bill::solvers::bsat2, false, true, false>; // last is false

  struct rewire_params : dependency::default_rewire_params
  {
    static constexpr uint32_t max_cuts_size = Params::max_cuts_size;
  };

  struct struct_params : dependency::default_struct_params
  {
    static constexpr uint32_t num_vars_sign = Params::max_num_leaves;
    static constexpr uint32_t max_cuts_size = Params::max_cuts_size;
  };

  struct window_params : dependency::default_window_params
  {
    static constexpr uint32_t max_num_leaves = Params::max_num_leaves;
    static constexpr uint32_t max_cuts_size = Params::max_cuts_size;
    static constexpr uint32_t max_cube_spfd = Params::max_cube_spfd;
  };

  using rewire_dependencies_t = dependency::rewire_dependencies<Ntk, rewire_params>;
  using struct_dependencies_t = dependency::struct_dependencies<Ntk, struct_params>;
  using window_dependencies_t = dependency::window_dependencies<Ntk, window_params>;
  // using simula_dependency_t = simula_dependency<Ntk, custom_simula_params>;

public:

  resynthesize_impl( Ntk& ntk, Database& database, Params ps, resynthesis_stats& st)
    : ntk_( ntk ),
      win_manager_( ntk, ps.window_manager_ps, st.window_st ),
      profiler_( ntk, win_manager_, ps.profiler_ps ),
      database_( database ),
      win_simulator_( ntk ),
      chain_simulator_( database.get_library() ),
      ps_( ps ),
      st_( st )
  {}
public:
  void run()
  {
    rewire_dependencies_t rewire_dependencies( ntk_ );
    struct_dependencies_t struct_dependencies( ntk_ );
    window_dependencies_t window_dependencies( ntk_ );

    chain_t best_chain;
    best_chain.add_inputs( Params::max_cuts_size );

    auto nmax = ntk_.size();
    profiler_.foreach_gate( [&]( auto n ) {
      if ( st_.num_optims >= ps_.stop_after_nth )
        return false;

      if ( n >= nmax )
        return false;
      /* Skip nodes which cannot result in optimization */
      if ( skip_node( n ) )
        return true;

      /* Build and run analysis on window */
      window_analysis( n );
      profiler_.init();

      if ( !win_manager_.is_valid() )
        return true;

      if ( ps_.try_rewire )
      {
        rewire_dependencies.run( win_manager_, win_simulator_ );
        std::optional<cut_t> best_cut;
        double best_reward = 0;

        rewire_dependencies.foreach_cut( [&]( auto& cut, auto i ) {
          auto const& cut_leaves = cut.leaves;
          auto const reward = profiler_.evaluate_rewiring( n, cut_leaves );
          std::cout << n << " rewire " << reward << " | best : " << best_reward << std::endl;
          if ( reward > best_reward )
          {
            best_reward = reward;
            best_cut = std::make_optional( cut );
          }
        } );
        if ( best_cut )
        {
          st_.num_optims++;
          auto const ids = ntk_.get_binding_ids( n );
          auto const fnew = ntk_.template create_node<Params::do_strashing>( ( *best_cut ).leaves, ids );
          substitute_node( n, fnew );
          return true;
        }
      }

      double best_reward = 0;
      std::vector<signal_t> best_leaves;

      if ( ps_.try_struct )
      {
        struct_dependencies.run( win_manager_, win_simulator_ );

        struct_dependencies.foreach_cut( [&]( auto& cut, auto i ) {
          best_reward = std::max( evaluate( cut, best_chain, best_leaves ), best_reward );
          std::cout << n << " struct best : " << best_reward << std::endl;
        } );

        if ( best_reward > 0 )
        {
          st_.num_optims++;
          auto const fnew = insert( ntk_, best_leaves, best_chain );
          substitute_node( n, fnew );
          return true;
        }
      }
      if ( ps_.try_window )
      {
        window_dependencies.run( win_manager_, win_simulator_ );
        window_dependencies.foreach_cut( [&]( auto& cut, auto i ) {
          best_reward = std::max( evaluate( cut, best_chain, best_leaves ), best_reward );
          std::cout << n << " window best : " << best_reward << std::endl;
        } );

        if ( best_reward > 0 )
        {
          st_.num_optims++;
          auto const fnew = insert( ntk_, best_leaves, best_chain );
          substitute_node( n, fnew );
          return true;
        }
      }

      return true;
    } );
  }

private:
  double evaluate( cut_t const& cut, chain_t& best_chain, std::vector<signal_t>& best_leaves )
  {
    double best_reward = 0;
    if ( cut.func.size() > 1 )
      return 0;
    double const cost_curr = profiler_.evaluate( cut.root, cut.leaves, cut.root );

    auto const cut_func = cut.func[0];
    auto signals = cut.leaves;
    auto times = get_times( signals );
    signals.resize( Params::max_cuts_size, std::numeric_limits<uint64_t>::max() );
    times.resize( Params::max_cuts_size, std::numeric_limits<double>::max() );
    signal_t best_signal;
    if ( !decomposer_.run( cut_func, times ) )
      return 0;

    bool const success = decomposer_.foreach_spec( [&]( auto& specs, uint8_t lit ) {
      auto const res = local_synthesis( cut, specs, lit, signals, times );
      if ( res )
      {
        auto const& [signal, time, sim] = *res;
        signals.push_back( signal );
        times.push_back( time );
        specs[lit].sim._bits = sim;
        best_signal = signal;
        return true;
      }
      return false;
    } );

    if ( !success )
      return 0;

    chain_t new_chain;
    new_chain.add_inputs( cut.leaves.size() );
    extract( new_chain, ntk_, cut.leaves, best_signal );
    double const cost_cand = profiler_.evaluate( new_chain, cut.leaves, cut.root );
    auto const reward = cost_curr - cost_cand;
    if ( reward > best_reward )
    {
      best_reward = reward;
      best_leaves = cut.leaves;
      best_chain = new_chain;
    }
    return best_reward;
  }

  std::optional<std::tuple<signal_t, double, func_t>> local_synthesis( cut_t const& cut, typename decomposer_t::specs_t const& specs, uint8_t lit, std::vector<signal_t>& signals, std::vector<double>& times )
  {
    std::vector<func_t const*> sim_ptrs;
    auto& spec = specs[lit];
    for ( auto i : spec.inputs )
      sim_ptrs.push_back( &specs[i].sim._bits );

    auto itt = dependency::extract_function<func_t, Database::max_num_vars>( sim_ptrs, specs[lit].sim._bits, specs[lit].sim._care );
    double best_loc_cost = std::numeric_limits<double>::max();
    std::optional<node> best_database_node;
    std::vector<signal> best_loc_leaves;
    std::vector<typename decomposer_t::cut_func_t const*> best_loc_sims;
    signal_t best_signal;

    enumerator_.foreach_dont_care_assignment( itt, sim_ptrs.size(), [&]( auto const& ctt ) {
      std::vector<signal> loc_leaves( spec.inputs.size() );
      std::vector<double> loc_times( spec.inputs.size() );
      std::vector<typename decomposer_t::cut_func_t const*> loc_sims = sim_ptrs;
      std::transform( spec.inputs.begin(), spec.inputs.end(), loc_leaves.begin(), [&]( auto const& lit ) { return signals[lit]; } );
      std::transform( spec.inputs.begin(), spec.inputs.end(), loc_times.begin(), [&]( auto const& lit ) { return times[lit]; } );
      // Perform boolean matching
      auto row = database_.boolean_matching( ctt, loc_times, loc_leaves, loc_sims );

      if ( row )
      {
        auto [index, cost_cand] = evaluate( *row, loc_leaves );
        if ( cost_cand < best_loc_cost )
        {
          best_loc_cost = cost_cand;
          best_database_node = std::make_optional( index );
          best_loc_leaves = loc_leaves;
          best_loc_sims = loc_sims;
        }
      }
    } );
    if ( best_database_node )
    {
      auto nnew = database_.write( *best_database_node, ntk_, best_loc_leaves );
      best_signal = ntk_.make_signal( nnew );
      signals.push_back( best_signal );
      times.push_back( profiler_.get_arrival( best_signal ) );
      chain_t loc_chain;
      loc_chain.add_inputs( Database::max_num_vars );
      extract( loc_chain, ntk_, best_loc_leaves, signals.back() );
      chain_simulator_( loc_chain, best_loc_sims );
      auto const sim = chain_simulator_.get_simulation( loc_chain, best_loc_sims, loc_chain.po_at( 0 ) );
      return std::make_optional( std::make_tuple( signals.back(), times.back(), sim ) );
    }
    return std::nullopt;
  }

  std::tuple<node_index_t, double> evaluate( uint64_t const& row, std::vector<signal_t> const& loc_leaves )
  {
    node_index_t best_database_node = std::numeric_limits<node_index_t>::max();
    double best_cost = std::numeric_limits<double>::max();

    database_.foreach_entry( row, [&]( auto const& entry ) {
      auto nnew = database_.write( entry, ntk_, loc_leaves );
      auto cost_new = profiler_.evaluate( nnew, loc_leaves, win_manager_.get_pivot() );
      if ( cost_new < best_cost )
      {
        best_cost = cost_new;
        best_database_node = entry.index;
      }
      ntk_.take_out_node( nnew );
    } );
    return std::make_tuple( best_database_node, best_cost );
  }

  void substitute_node( node_t const& n, signal_t const& fnew )
  {
    auto const nnew = ntk_.get_node( fnew );
    std::vector<signal_t> fs;
    ntk_.foreach_output( nnew, [&]( auto f ) {
      fs.push_back( f );
    } );
    ntk_.substitute_node( n, fs );
  }

private:
  /*! \brief Checks if the node should be analyzed for optimization or skipped */
  bool skip_node( node const& n )
  {
    if ( ntk_.fanout_size( n ) > ps_.fanout_limit )
      return true;
    if ( ntk_.fanout_size( n ) <= 0 )
      return true;
    if ( ntk_.is_pi( n ) || ntk_.is_constant( n ) )
      return true;
    if ( ntk_.is_dead( n ) )
      return true;
    return false;
  }

  /*! \brief Perform window analysis if required by any heuristic. Return false if failure */
  void window_analysis( node const& n )
  {
    win_manager_.run( n );
    win_simulator_.run( win_manager_ );
  }

  std::vector<double> get_times( std::vector<signal> const& leaves )
  {
    assert( Profiler::has_arrival && "[e] The profiler does not have the arrival tracker" );
    std::vector<double> times( leaves.size() );
    std::transform( leaves.begin(), leaves.end(), times.begin(),
                    [&]( auto const& f ) { return profiler_.get_arrival( f ); } );
    return times;
  }

private:
  Ntk& ntk_;
  dependency::function_enumerator<Database::max_num_vars> enumerator_;
  window_manager_t win_manager_;
  windowing::window_simulator<Ntk, Params::window_manager_params::max_num_leaves> win_simulator_;
  Profiler profiler_;
  Database& database_;
  decomposer_t decomposer_;
  evaluation::chain_simulator<chain_t, func_t> chain_simulator_;
  Params ps_;
  resynthesis_stats& st_;
};

} /* namespace detail */

template<class Ntk, class Database, typename Params = default_resynthesis_params<RINOX_MAX_NUM_LEAVES>>
void area_resynthesize( Ntk& ntk, Database& database, Params ps = {}, resynthesis_stats* pst = nullptr )
{
  using WinMngr = windowing::window_manager<Ntk, typename Params::window_manager_params>;
  using Profiler = profilers::area_profiler<Ntk, WinMngr>;
  resynthesis_stats st;
  detail::resynthesize_impl<Ntk, Database, Profiler, Params> p( ntk, database, ps, st );
  p.run();
  if ( pst != nullptr )
    *pst = st;
}

template<class Ntk, class Database, typename Params = default_resynthesis_params<RINOX_MAX_NUM_LEAVES>>
void delay_resynthesize( Ntk& ntk, Database& database, Params ps = {}, resynthesis_stats* pst = nullptr )
{
  using WinMngr = windowing::window_manager<Ntk, typename Params::window_manager_params>;
  using Profiler = profilers::delay_profiler<Ntk, WinMngr>;
  resynthesis_stats st;
  detail::resynthesize_impl<Ntk, Database, Profiler, Params> p( ntk, database, ps, st );
  p.run();
  if ( pst != nullptr )
    *pst = st;
}

template<class Ntk, class Database, typename Params = default_resynthesis_params<RINOX_MAX_NUM_LEAVES>>
void power_resynthesize( Ntk& ntk, Database& database, Params ps = {}, resynthesis_stats* pst = nullptr )
{
  using WinMngr = windowing::window_manager<Ntk, typename Params::window_manager_params>;
  using Profiler = profilers::power_profiler<Ntk, WinMngr, Params::window_manager_params::max_num_leaves>;
  resynthesis_stats st;
  detail::resynthesize_impl<Ntk, Database, Profiler, Params> p( ntk, database, ps, st );
  p.run();
  if ( pst != nullptr )
    *pst = st;
}

} // namespace algorithms
} // namespace opto
} // namespace rinox