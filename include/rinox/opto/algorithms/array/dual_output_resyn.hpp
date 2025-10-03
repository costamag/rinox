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
  \file multiple_output.hpp
  \brief Resynthesis using multiple output gates

  \author Andrea Costamagna
*/

#pragma once

#include "../../../network/visualize.hpp"
#include "../../../windowing/window_manager.hpp"
#include "../../../windowing/window_simulator.hpp"
#include "../../../boolean/spfd.hpp"
#include "../../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../profilers/area_profiler.hpp"
#include "../../profilers/delay_profiler.hpp"
#include "../../profilers/power_profiler.hpp"
#include "../../profilers/profilers_utils.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <chrono>
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

#define DEBUG_INFO true

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
  double time{ 0 };

  /*! \brief Expected gain. */
  uint32_t estimated_gain{ 0 };

  /*! \brief Candidates */
  uint32_t candidates{ 0 };

  uint32_t num_struct{ 0 };
  uint32_t num_window{ 0 };
  uint32_t num_simula{ 0 };
  uint32_t num_rewire{ 0 };
  uint32_t num_optims{ 0 };
  double fDual{ 0 };

  void report() const
  {
    std::cout << fmt::format( "    num struct       = {:5d}\n", num_struct );
    std::cout << fmt::format( "    num window       = {:5d}\n", num_window );
    std::cout << fmt::format( "    num simula       = {:5d}\n", num_simula );
    std::cout << fmt::format( "    num rewire       = {:5d}\n", num_rewire );
  }
};

template<uint32_t NumLeaves = RINOX_MAX_NUM_LEAVES>
struct default_dual_output_resyn_params
{
  static constexpr uint32_t max_num_leaves = NumLeaves;

  profilers::profiler_params profiler_ps;

  /*! \brief If true, candidates are only accepted if they do not increase logic depth. */
  struct window_manager_params : windowing::default_window_manager_params
  {
    static constexpr uint32_t max_num_leaves = NumLeaves;
    bool preserve_depth = true;
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
  static constexpr uint32_t max_num_vars_dual = 5u;
  static constexpr uint32_t max_num_vars_single = 6u;

  /*! \brief Maximum fanout size for a node to be optimized*/
  uint32_t fanout_limit = 12u;
};

namespace detail
{

template<class Ntk, typename Profiler, typename Params = default_dual_output_resyn_params<RINOX_MAX_NUM_LEAVES>>
class area_dual_output_resyn_impl
{
  using node = typename Ntk::node;
  using node_t = typename Ntk::node;
  using signal = typename Ntk::signal;
  using signal_t = typename Ntk::signal;
  using node_index_t = typename Ntk::node;
  using window_manager_t = windowing::window_manager<Ntk, typename Params::window_manager_params>;
  using signals_t = rinox::math::static_vector<signal_t, Params::max_cuts_size>;
  using times_t = rinox::math::static_vector<double, Params::max_cuts_size>;
  static constexpr uint32_t max_num_leaves = Params::window_manager_params::max_num_leaves;
  static constexpr uint32_t K1 = Params::max_num_vars_single;
  static constexpr uint32_t K2 = Params::max_num_vars_dual;
  using signature_t = kitty::static_truth_table<max_num_leaves>;
  using clock = std::chrono::high_resolution_clock;

public:
  area_dual_output_resyn_impl( Ntk& ntk, Params ps, resynthesis_stats& st )
      : ntk_( ntk ),
        arrival_( ntk_ ),
        win_manager_( ntk, ps.window_manager_ps, st.window_st ),
        win_simulator_( ntk ),
        ps_( ps ),
        st_( st )
  {
    for ( int n = 0; n < 7; ++n )
    {
      projs.emplace_back();
      for ( int i = 0; i < n; ++i ){
        projs[n].emplace_back( n );
        kitty::create_nth_var( projs[n][i], i );
      }
    }
  }

public:
  void run()
  {
    auto start = clock::now();
    try_merge();
    //try_single_0_resub();
    //try_dual_0_resub();
    //try_single_1_resub();
    //try_dual_1_resub();
    auto stop = clock::now();
    st_.time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count());

    double total = 0;
    st_.fDual = 0;
    ntk_.foreach_gate( [&]( auto const& n ){
      if ( ntk_.num_outputs( n ) == 2 ) st_.fDual+=1;
      total+=1;
    } );
    st_.fDual/=total;
  }

private:

void try_dual_1_resub()
{
  const auto nmax = ntk_.size();
  int cnt = 0;

  ntk_.foreach_gate([&](auto const& n) {
    if ( ntk_.is_pi(n) ) return true;
    if (cnt++ >= nmax) return false;
    if (st_.num_optims >= ps_.stop_after_nth) return false;
    if (skip_node(n)) return true;
    if (ntk_.num_outputs(n) != 2) return true;

    win_manager_.run( n );
    if ( win_manager_.num_mffc() <= 1 ) return true;
    win_simulator_.run( win_manager_ );

    // check if a divisor has the same functionality
    std::vector<signature_t const *> sims;
    ntk_.foreach_output( n, [&]( auto const& f ){ sims.push_back( &win_simulator_.get(f) ); } );
    spfds_.init( sims );
    std::vector<signal> children;
    if ( spfds_.is_covered() ) return true;
    while ( children.size() < K2 && !spfds_.is_covered() )
    {
      std::optional<signal> best_div;
      uint32_t best_cost = spfds_.get_num_edges();
      uint32_t cost;
      win_manager_.foreach_divisor( [&]( auto d, auto i ){
        auto const nd = ntk_.get_node( d );
        cost = spfds_.evaluate( win_simulator_.get( d ) );
        if ( cost < best_cost )
        {
          best_cost = cost;
          best_div = d;
        }
      } );
      if ( best_div )
      {
        children.push_back( *best_div );
        spfds_.update( win_simulator_.get( *best_div ) );
      }
    }
    if ( ps_.window_manager_ps.preserve_depth )
    {
      double time = 0;
      for ( auto const f : children )
        time = std::max( time, arrival_.get_time( f ) );
      if ( (time + 1) > arrival_.get_time( ntk_.make_signal( n ) ) ) return true;
    }
    if ( spfds_.is_covered() )
    {
      auto tt0 = extract_function( children, signal{ n, 0 } );
      auto tt1 = extract_function( children, signal{ n, 1 } );
      std::vector<kitty::dynamic_truth_table> tts = { tt0, tt1 };
      auto fnew = ntk_.create_node( children, tts );
      std::vector<signal> new_signals { signal{ fnew.index, 0 }, signal{ fnew.index, 1 } };
      ntk_.substitute_node( n, new_signals );
      return true;
    }

    return true;
  });
}

void try_single_1_resub()
{
  const auto nmax = ntk_.size();
  int cnt = 0;

  ntk_.foreach_gate([&](auto const& n) {
    if ( ntk_.is_pi(n) ) return true;
    if (cnt++ >= nmax) return false;
    if (st_.num_optims >= ps_.stop_after_nth) return false;
    if (skip_node(n)) return true;
    if (ntk_.num_outputs(n) != 1) return true;

    win_manager_.run( n );
    if ( win_manager_.num_mffc() <= 1 ) return true;
    win_simulator_.run( win_manager_ );

    // check if a divisor has the same functionality
    std::vector<signature_t const *> sims;
    ntk_.foreach_output( n, [&]( auto const& f ){ sims.push_back( &win_simulator_.get(f) ); } );
    spfds_.init( sims );

    std::optional<signal> best_div;
    uint32_t best_cost = std::numeric_limits<uint32_t>::max();
    uint32_t cost;
    win_manager_.foreach_divisor( [&]( auto d, auto i ){
      auto const nd = ntk_.get_node( d );
      cost = spfds_.evaluate( win_simulator_.get( d ) );
      if ( cost < best_cost )
      {
        best_cost = cost;
        best_div = d;
      }
    } );

    if ( win_manager_.num_mffc() > 1 )
    {
      std::vector<signal> children = { *best_div };

      while ( children.size() < K1 && !spfds_.is_covered() )
      {
        std::optional<signal> best_div;
        uint32_t best_cost = spfds_.get_num_edges();
        uint32_t cost;
        win_manager_.foreach_divisor( [&]( auto d, auto i ){
          auto const nd = ntk_.get_node( d );
          cost = spfds_.evaluate( win_simulator_.get( d ) );
          if ( cost < best_cost )
          {
            best_cost = cost;
            best_div = d;
          }
        } );
        if ( best_div )
        {
          children.push_back( *best_div );
          spfds_.update( win_simulator_.get( *best_div ) );
        }
      }
      if ( ps_.window_manager_ps.preserve_depth )
      {
        double time = 0;
        for ( auto const f : children )
          time = std::max( time, arrival_.get_time( f ) );
        if ( (time + 1) > arrival_.get_time( ntk_.make_signal( n ) ) ) return true;
      }
      if ( spfds_.is_covered() )
      {
        auto tt = extract_function( children, ntk_.make_signal( n ) );
        auto fnew = ntk_.create_node( children, tt );

        ntk_.substitute_node( n, fnew );
        return true;
      }
    }

    return true;
  });
}

void try_dual_0_resub()
{
  const auto nmax = ntk_.size();
  int cnt = 0;

  ntk_.foreach_gate([&](auto const& n) {
    if ( ntk_.is_pi(n) ) return true;
    if (cnt++ >= nmax) return false;
    if (st_.num_optims >= ps_.stop_after_nth) return false;
    if (skip_node(n)) return true;
    if (ntk_.num_outputs(n) != 1) return true;

    window_analysis(n);
    // check if a divisor has the same functionality
    std::vector<signature_t const *> sims;
    ntk_.foreach_output( n, [&]( auto const& f ){ sims.push_back( &win_simulator_.get(f) ); } );
    spfds_.init( sims );

    std::optional<signal> best_div;
    uint32_t best_cost = std::numeric_limits<uint32_t>::max();
    uint32_t cost;
    win_manager_.foreach_divisor( [&]( auto d, auto i ){
      auto const nd = ntk_.get_node( d );
      cost = spfds_.evaluate( win_simulator_.get( d ) );
      if ( cost < best_cost )
      {
        best_cost = cost;
        best_div = d;
      }
    } );

     // try constructing an ad-hoc dual-output node
    if ( best_div &&
         ( !win_manager_.is_input( ntk_.get_node( *best_div ) ) ) &&
         ( ntk_.fanin_size( *best_div ) < K2 ) &&
         ( ntk_.get_node( *best_div ) == 1 )
      )
    {
      auto const nd = ntk_.get_node( *best_div );
      ntk_.foreach_fanin( nd, [&]( auto const& fi, auto ii ){
        spfds_.update( win_simulator_.get( fi ) );
      } );
      std::vector<signal> children = ntk_.get_children( nd );
      while ( children.size() < K2 && !spfds_.is_covered() )
      {
        std::optional<signal> best_div;
        uint32_t best_cost = std::numeric_limits<uint32_t>::max();
        uint32_t cost;
        win_manager_.foreach_divisor( [&]( auto d, auto i ){
          auto const nd = ntk_.get_node( d );
          cost = spfds_.evaluate( win_simulator_.get( d ) );
          if ( cost < best_cost )
          {
            best_cost = cost;
            best_div = d;
          }
        } );
        if ( best_div )
        {
          children.push_back( *best_div );
          spfds_.update( win_simulator_.get( *best_div ) );
        }
      }

      if ( spfds_.is_covered() )
      {
        if ( ps_.window_manager_ps.preserve_depth )
        {
          double time = 0;
          for ( auto const f : children )
            time = std::max( time, arrival_.get_time( f ) );
          if ( (time + 1) > arrival_.get_time( ntk_.make_signal( n ) ) ) return true;
        }
        auto tt = extract_function( children, ntk_.make_signal( n ) );
        auto fnew = ntk_.create_node( children, tt );
        merge( fnew.index, nd );
        return true;
      }
    }
    return true;
  });
}


kitty::dynamic_truth_table extract_function( std::vector<signal> const& children, signal f )
{
  signature_t func = win_simulator_.get( f );
  kitty::dynamic_truth_table tt = projs[children.size()][0].construct();
  signature_t tmp_s;
  kitty::dynamic_truth_table tmp_f( children.size() );
  int num_minterms = 1u << children.size();
  for ( int m = 0; m < num_minterms; ++m )
  {
    tmp_s |= ~tmp_s;
    for ( int i = 0; i < children.size(); ++i )
    {
      if ( ( ( m >> i ) & 1u ) > 0 )
      {
        tmp_s &= win_simulator_.get( children[i] );
      }
      else
      {
        tmp_s &= win_simulator_.get( children[i] );
      }
    }
    if ( kitty::count_ones( tmp_s & func ) > 0 )
      kitty::set_bit( tt, m );
    else
      kitty::clear_bit( tt, m );
  }
  return tt;
}

// either a node has the same signature as a divisor
void try_single_0_resub()
{
  const auto nmax = ntk_.size();
  int cnt = 0;

  ntk_.foreach_gate([&](auto const& n) {
    if ( ntk_.is_pi(n) ) return true;
    if (cnt++ >= nmax) return false;
    if (st_.num_optims >= ps_.stop_after_nth) return false;
    if (skip_node(n)) return true;
    if (ntk_.num_outputs(n) != 1) return true;

    window_analysis(n);
    // check if a divisor has the same functionality
    signal f = ntk_.make_signal( n );
    auto const& sim_root = win_simulator_.get( f );
    std::optional<std::pair<signal, bool>> best;
    win_manager_.foreach_divisor( [&]( auto const& d, auto i ){
      if ( !ps_.window_manager_ps.preserve_depth || ( arrival_.get_time( f ) >= arrival_.get_time( d ) ) )
      {
        auto const& sim_div = win_simulator_.get( d );
        if ( kitty::equal( sim_div, sim_root ) )
        {
          best = std::pair{d, false};
        }
        else if ( kitty::equal( ~sim_div, sim_root ) && !best )
        {
          best = std::pair{d, true};
        }
      }
    } );
    if ( best )
    {
      if ( (*best).second )
      {
        kitty::dynamic_truth_table inv(1u);
        kitty::create_from_binary_string( inv, "01" );
        auto fnew = ntk_.create_node( { (*best).first }, inv );
        ntk_.substitute_node( n, fnew );
      }
      else
      {
        ntk_.substitute_node( n, (*best).first );
      }
      st_.num_optims++;
    }
    return true;
  });
}

bool reachable(node_index_t src, node_index_t dst)
{
  auto fsrc = ntk_.make_signal( src );
  auto fdst = ntk_.make_signal( dst );
  double min_time = std::min( arrival_.get_time( fsrc ), arrival_.get_time( fdst ) )-2;
  // Is there a path from src to dst?
  ntk_.incr_trav_id();
  bool found = false;
  std::function<void(node)> dfs = [&](node u) {
    if (found) return;
    if (u == dst) { found = true; return; }
    ntk_.foreach_fanin(u, [&](auto fi, auto) {
      auto v = ntk_.get_node(fi);
      if (ntk_.visited(v) == ntk_.trav_id() || ( arrival_.get_time( fi ) < min_time )) return;
      ntk_.set_visited(v, ntk_.trav_id());
      dfs(v);
    });
  };
  dfs(src);
  return found;
};

void try_merge()
{
  const auto nmax = ntk_.size();
  int cnt = 0;

  ntk_.foreach_gate([&](auto const& n) {
    if ( ntk_.is_pi(n) ) return true;
    if (cnt++ >= nmax) return false;
    if (st_.num_optims >= ps_.stop_after_nth) return false;
    if (skip_node(n)) return true;

    if (ntk_.num_outputs(n) != 1) return true;

    // Mark fanouts of n
    const auto color_fo = (n << 2u) | 2u;
    ntk_.foreach_fanout(n, [&](auto const& no) {
      ntk_.set_visited(no, color_fo);
    });

    // Collect inputs of n
    std::set<signal_t> base_inputs;
    ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
      base_inputs.insert(fi);
    });

    std::optional<node_index_t> best_node;
    int best_cost = std::numeric_limits<int>::max();

    // Explore fanouts of fanins
    ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
      ntk_.foreach_fanout(fi.index, [&](auto const& no) {
        if (n == no) return;
        if (ntk_.num_outputs(no) != 1u) return;
        if (ntk_.visited(no) == color_fo) return;

        // Compute combined inputs size
        std::set<signal_t> merged_inputs = base_inputs;
        ntk_.foreach_fanin(no, [&](auto const& foi, auto) {
          merged_inputs.insert(foi);
        });

        if ( !reachable(n, no) && !reachable(no, n) && ( (int)merged_inputs.size() < best_cost)) {
          if ( !ps_.window_manager_ps.preserve_depth || ( std::abs(arrival_.get_time( ntk_.make_signal( n ) ) - arrival_.get_time( ntk_.make_signal( no ) )) < 1 ) )
          {
            best_node = no;
            best_cost = (int)merged_inputs.size();
          }
        }
      });
    });

    if (best_node) {
      merge(n, *best_node);
    }

    return true;
  });
}

  void merge( node_index_t n1, node_index_t n2 )
  {
    if ( ntk_.is_pi( n1 ) || ntk_.is_pi( n2 ) ) return;
    auto tt1 = ntk_.node_function( n1 );
    auto tt2 = ntk_.node_function( n2 );

    std::unordered_map<uint64_t, std::pair<int, int>> map;
    std::set<signal_t> inps;
    ntk_.foreach_fanin( n1, [&]( auto const& fi, auto ii ){ inps.insert( fi ); } );
    ntk_.foreach_fanin( n2, [&]( auto const& fi, auto ii ){ inps.insert( fi ); } );
    if ( inps.size() > K2 ) return;
    int i = 0;
    std::vector<signal_t> children;
    for ( auto const& f : inps )
    {
      children.push_back( f );
      map[f.data] = std::pair{ inps.size(), i++ };
    }
    std::vector<kitty::dynamic_truth_table const*> sim_ptrs1;
    std::vector<kitty::dynamic_truth_table const*> sim_ptrs2;

    ntk_.foreach_fanin( n1, [&]( auto const& fi, auto ii ){
      auto const p = map[fi.data];
      sim_ptrs1.push_back( &projs[p.first][p.second] );
    } );
    ntk_.foreach_fanin( n2, [&]( auto const& fi, auto ii ){
      auto const p = map[fi.data];
      sim_ptrs2.push_back( &projs[p.first][p.second] );
    } );
    auto funcs1 = ntk_.compute( n1, sim_ptrs1 );
    auto funcs2 = ntk_.compute( n2, sim_ptrs2 );
    std::vector<kitty::dynamic_truth_table> tts = { funcs1[0], funcs2[0] };

    auto fnew = ntk_.create_node( children, tts );
    ntk_.substitute_node( n1, signal_t{ fnew.index, 0 } );
    ntk_.substitute_node( n2, signal_t{ fnew.index, 1 } );
    st_.num_optims++;
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


private:
  Ntk& ntk_;
  window_manager_t win_manager_;
  windowing::window_simulator<Ntk, Params::window_manager_params::max_num_leaves> win_simulator_;
  Params ps_;
  analyzers::trackers::arrival_times_tracker<Ntk> arrival_;
  resynthesis_stats& st_;
  boolean::spfds_manager<signature_t, (1u << K2), 2u> spfds_;
  std::vector<std::vector<kitty::dynamic_truth_table>> projs;

};

} /* namespace detail */

template<class Ntk, typename Params = default_dual_output_resyn_params<RINOX_MAX_NUM_LEAVES>>
void area_dual_output_resyn( Ntk& ntk, Params ps = {}, resynthesis_stats* pst = nullptr )
{
  using WinMngr = windowing::window_manager<Ntk, typename Params::window_manager_params>;
  using Profiler = profilers::area_profiler<Ntk, WinMngr>;
  resynthesis_stats st;
  detail::area_dual_output_resyn_impl<Ntk, Profiler, Params> p( ntk, ps, st );

  p.run();
  if ( pst != nullptr )
    *pst = st;
}

} // namespace algorithms
} // namespace opto
} // namespace rinox