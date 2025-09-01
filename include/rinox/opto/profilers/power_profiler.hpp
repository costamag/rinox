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
  \file power_profiler.hpp
  \brief Analyzer of the power

  \author Andrea Costamagna
*/

#pragma once

#include "../../analyzers/analyzers_utils/switching.hpp"
#include "../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../analyzers/trackers/gate_load_tracker.hpp"
#include "../../analyzers/trackers/sensing_times_tracker.hpp"
#include "../../databases/mapped_database.hpp"
#include "profilers_utils.hpp"

namespace rinox
{

namespace opto
{

namespace profilers
{

template<class Ntk, typename WinMngr, uint32_t MaxNumLeaves = 12u>
class power_profiler
{
public:
  using func_t = kitty::static_truth_table<MaxNumLeaves>;
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;
  using cost_t = double;
  static cost_t constexpr min_cost = std::numeric_limits<cost_t>::min();
  static cost_t constexpr max_cost = std::numeric_limits<cost_t>::max();
  static bool constexpr pass_window = true;
  static bool constexpr has_arrival = true;
  static constexpr uint32_t max_num_steps = 10u;

  struct node_with_cost_t
  {
    node_index_t root;
    cost_t mffc_power;
  };

public:
  power_profiler( Ntk& ntk, WinMngr& win_manager, profiler_params const& ps )
      : ntk_( ntk ),
        ps_( ps ),
        nodes_( ntk_.size() ),
        workload_( MaxNumLeaves ),
        arrival_( ntk ),
        loading_( ntk ),
        sensing_( ntk ),
        signal_to_activity_( ntk ),
        win_manager_( win_manager )
  {
  }

  void init()
  {
    // assign the simulation to the input signals
    assert( win_manager_.is_valid() );
    signal_to_activity_.resize();
    activity_.reserve( win_manager_.size() );
    activity_.resize( MaxNumLeaves );
    win_manager_.foreach_input( [&]( auto const& f, auto const& i ) {
      assert( i < MaxNumLeaves );
      signal_to_activity_[f] = i;
      activity_[i] = workload_.get( i );
    } );
    simulate_window();
  }

  void simulate_window()
  {

    /* simulate the network in topological order */
    win_manager_.foreach_divisor( [&]( auto f, auto i ) {
      (void)i;
      if ( !win_manager_.is_input( ntk_.get_node( f ) ) )
        simulate_signal( f );
    } );
    win_manager_.foreach_mffc( [&]( auto n, auto i ) {
      (void)i;
      ntk_.foreach_output( n, [&]( auto const& f ) {
        simulate_signal( f );
      } );
    } );
    win_manager_.foreach_tfo( [&]( auto n, auto i ) {
      (void)i;
      ntk_.foreach_output( n, [&]( auto const& f ) {
        simulate_signal( f );
      } );
    } );
  }

  void simulate_signal( signal_t const& f )
  {
    auto const& children = ntk_.get_children( f );
    simulate_signal( f, children );
  }

  template<typename Signals>
  void simulate_signal( signal_t const& f, Signals const& fanin )
  {
    double const norm = static_cast<double>( workload_.num_bits() );

    std::vector<func_t const*> sim_ptrs;
    auto const binding = ntk_.get_binding( f );
    /* simulate the beginning of the clock cycle */
    signal_to_activity_[f] = activity_.size();
    activity_.emplace_back();
    sim_ptrs.clear();
    for ( auto const& fi : fanin )
    {
      sim_ptrs.push_back( &activity_[signal_to_activity_[fi]][0u] );
    }
    ntk_.compute( activity_.back()[0u], f, sim_ptrs );

    uint32_t step = 1;
    for ( ; step < max_num_steps; ++step )
      activity_.back()[step] = activity_.back()[0u];

    /* simulate the end of the clock cycle */
    sim_ptrs.clear();
    for ( auto const& fi : fanin )
    {
      auto const index = signal_to_activity_[fi];
      sim_ptrs.push_back( &activity_[index][max_num_steps - 1] );
    }
    ntk_.compute( activity_.back()[max_num_steps - 1], f, sim_ptrs );

    step = 1;
    if ( arrival_.get_time( f ) > sensing_.get_time( f ) )
    {
      while ( step < ( max_num_steps - 1 ) )
      {
        double const time = get_time( step, sensing_.get_time( f ) - binding.avg_pin_delay, arrival_.get_time( f ) + binding.avg_pin_delay );
        sim_ptrs.clear();
        for ( auto ii = 0u; ii < fanin.size(); ++ii )
        {
          auto const fi = fanin[ii];
          assert( binding.max_pin_time.size() >= fanin.size() );

          double const time_i = time - binding.max_pin_time[ii];
          auto step_i = get_step( time_i, sensing_.get_time( fi ) - binding.avg_pin_delay, arrival_.get_time( fi ) + binding.avg_pin_delay );
          auto const index = signal_to_activity_[fi];
          sim_ptrs.push_back( &activity_[index][step_i] );
        }
        ntk_.compute( activity_.back()[step], f, sim_ptrs );

        step++;
      }
    }
    else
    {
      /* just replicate the simulation */
      while ( step++ < max_num_steps / 2 )
        activity_.back()[step] = activity_.back()[0u];
    }
    while ( step++ < ( max_num_steps - 1 ) )
    {
      activity_.back()[step] = activity_.back()[max_num_steps - 1];
    }
    double glitching = 0;
    double switching = 0;
    double zerodelay = 0;

    for ( auto step = 1u; step < max_num_steps; ++step )
      switching += kitty::count_ones( activity_.back()[step] ^ activity_.back()[step - 1] );

    zerodelay = kitty::count_ones( activity_.back()[0] ^ activity_.back()[max_num_steps - 1] );

    glitching = switching - zerodelay;
    /* normalize */
    glitching /= norm;
    switching /= norm;
    activity_.back().set_glitching( glitching );
    activity_.back().set_switching( switching );
    activity_.back().set_dyn_power( loading_.get_load( f ) * switching );
  }

  double get_arrival( signal_t const& f ) const
  {
    return arrival_.get_time( f );
  }

  template<class List_t, typename Signals>
  cost_t evaluate( List_t const& list, Signals const& leaves, node_index_t const nold )
  {
    auto const size_before = activity_.size();
    activity_.reserve( size_before + list.num_gates() );

    signal_t const f = insert( ntk_, leaves, list );
    signal_to_activity_.resize();

    simulate_tfi( f, leaves );

    node_index_t const n = ntk_.get_node( f );
    cost_t cost_deref = recursive_deref( n );
    cost_t cost_ref = recursive_ref( n );
    assert( std::abs( cost_ref - cost_deref ) < ps_.eps && "[e] Cost ref and deref should be the same" );

    ntk_.foreach_fanout( nold, [&]( auto const& no ) {
      ntk_.foreach_fanin( no, [&]( auto const& fi, auto ii ) {
        if ( ntk_.get_node( fi ) == nold )
        {
          auto const switching = activity_[signal_to_activity_[fi]].get_switching();
          auto const load = ntk_.get_input_load( f, ii );
          cost_deref += load * switching;
        }
      } );
    } );

    if ( ntk_.fanout_size( n ) == 0 )
      ntk_.take_out_node( n );

    activity_.resize( size_before );

    return cost_deref;
  }

  template<typename Signals>
  cost_t evaluate_rewiring( node_index_t const& n, Signals const& new_children )
  {
    double cost_curr = 0.0;
    ntk_.foreach_output( n, [&]( auto const& f ) {
      auto const switching = activity_[signal_to_activity_[f]].get_switching();
      auto const load = loading_.get_load( f );
      cost_curr += switching * load;
    } );

    ntk_.foreach_fanin( n, [&]( auto const& fi, auto const ii ) {
      auto const switching = activity_[signal_to_activity_[fi]].get_switching();
      auto const load = ntk_.get_input_load( ntk_.make_signal( n ), ii );
      cost_curr += switching * load;
    } );

    double cost_cand = 0.0;
    for ( auto i = 0u; i < new_children.size(); ++i )
    {
      auto const& fi = new_children[i];
      auto const switching = activity_[signal_to_activity_[fi]].get_switching();
      auto const load = ntk_.get_input_load( ntk_.make_signal( n ), i );
      cost_cand += switching * load;
    }

    ntk_.foreach_output( n, [&]( auto const& f ) {
      auto const tmp = activity_[signal_to_activity_[f]];

      simulate_signal( f, new_children );
      auto const switching = activity_[signal_to_activity_[f]].get_switching();
      auto const load = loading_.get_load( f );
      cost_cand += switching * load;
      activity_[signal_to_activity_[f]] = tmp; // restore the activity
    } );

    return cost_curr - cost_cand;
  }

  template<typename Signals>
  cost_t evaluate( node_index_t const& n, Signals const& children, node_index_t const& nold )
  {
    auto const size_before = activity_.size();
    signal_to_activity_.resize();
    simulate_tfi( ntk_.make_signal( n ), children );

    std::vector<node_index_t> leaves( children.size() );
    std::transform( children.begin(), children.end(), leaves.begin(), [&]( auto const f ) { return ntk_.get_node( f ); } );
    cost_t cost_deref = measure_mffc_deref( n, leaves );
    cost_t cost_ref = measure_mffc_ref( n, leaves );
    assert( std::abs( cost_ref - cost_deref ) < ps_.eps && "[e] Cost ref and deref should be the same" );

    ntk_.foreach_fanout( nold, [&]( auto const& no ) {
      ntk_.foreach_output( no, [&]( auto const& f ) {
        ntk_.foreach_fanin( no, [&]( auto const& fi, auto ii ) {
          if ( ntk_.get_node( fi ) == nold )
          {
            auto const switching = activity_[signal_to_activity_[fi]].get_switching();
            auto const load = ntk_.get_input_load( f, ii );
            cost_deref += load * switching;
          }
        } );
      } );
    } );
    activity_.resize( size_before );
    return cost_deref;
  }

  template<typename Signals>
  void simulate_tfi( signal_t const& f, Signals const& leaves )
  {
    ntk_.incr_trav_id();
    for ( auto const& l : leaves )
      ntk_.set_visited( ntk_.get_node( l ), ntk_.trav_id() );
    simulate_rec( f, ntk_.trav_id() );
  }

  void simulate_rec( signal_t const& f, uint32_t color )
  {
    auto const n = ntk_.get_node( f );
    if ( ( ntk_.visited( n ) == color ) || ( win_manager_.is_contained( n ) ) )
      return;
    ntk_.set_visited( n, color );

    ntk_.foreach_fanin( f, [&]( auto const& fi, auto ii ) {
      (void)ii;
      simulate_rec( fi, color );
    } );
    simulate_signal( f, ntk_.get_children( n ) );
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn )
  {
    ntk_.foreach_gate( [&]( auto const& n ) {
      fn( n );
    } );
  }

private:
  /*! \brief get the simulation time given the simulation step in the activity window */
  double get_time( uint32_t step, double s, double a ) const
  {
    if ( a <= s )
      return s; // or a; pick a convention
    step = std::min( max_num_steps - 1, step );
    return s + step * ( a - s ) / ( max_num_steps - 1 );
  }

  uint32_t get_step( double t, double s, double a ) const
  {
    if ( a <= s )
      return 0u;
    const double sf = ( max_num_steps - 1 ) * ( t - s ) / ( a - s );
    if ( !std::isfinite( sf ) )
      return 0u;
    auto step = static_cast<uint32_t>( std::lround( sf ) );
    return std::min( max_num_steps - 1, step );
  }

  cost_t recursive_deref( node_index_t const& n )
  {
    /* terminate? */
    if ( ntk_.is_constant( n ) || ntk_.is_pi( n ) )
      return 0.0;

    /* recursively collect nodes */
    cost_t power = 0;

    ntk_.foreach_fanin( n, [&]( auto const& fi, auto i ) {
      node_index_t const ni = ntk_.get_node( fi );
      auto const load = ntk_.get_input_load( ntk_.make_signal( n ), i );
      power += activity_[signal_to_activity_[fi]].get_switching() * load;
      if ( ntk_.decr_fanout_size( ni ) == 0 )
      {
        power += recursive_deref( ni );
      }
    } );
    return power;
  }

  cost_t recursive_ref( node_index_t const& n )
  {
    /* terminate? */
    if ( ntk_.is_constant( n ) || ntk_.is_pi( n ) )
      return 0.0;

    /* recursively collect nodes */
    cost_t power = 0;

    ntk_.foreach_fanin( n, [&]( auto const& fi, auto i ) {
      node_index_t const ni = ntk_.get_node( fi );
      auto const load = ntk_.get_input_load( ntk_.make_signal( n ), i );
      power += activity_[signal_to_activity_[fi]].get_switching() * load;
      if ( ntk_.incr_fanout_size( ni ) == 0 )
      {
        power += recursive_ref( ni );
      }
    } );
    return power;
  }

  cost_t measure_mffc_deref( node_index_t const& n, std::vector<node_index_t> const& leaves )
  {
    /* reference cut leaves */
    for ( auto l : leaves )
    {
      if ( l < std::numeric_limits<uint32_t>::max() )
        ntk_.incr_fanout_size( l );
    }

    cost_t mffc_cost = recursive_deref( n );

    /* dereference leaves */
    for ( auto l : leaves )
    {
      if ( l < std::numeric_limits<uint32_t>::max() )
        ntk_.decr_fanout_size( l );
    }

    return mffc_cost;
  }

  cost_t measure_mffc_ref( node_index_t const& n, std::vector<node_index_t> const& leaves )
  {
    /* reference cut leaves */
    for ( auto l : leaves )
    {
      if ( l < std::numeric_limits<uint32_t>::max() )
        ntk_.incr_fanout_size( l );
    }

    cost_t mffc_cost = recursive_ref( n );

    /* dereference leaves */
    for ( auto l : leaves )
    {
      if ( l < std::numeric_limits<uint32_t>::max() )
        ntk_.decr_fanout_size( l );
    }

    return mffc_cost;
  }

private:
  Ntk& ntk_;
  profiler_params const& ps_;
  std::vector<node_with_cost_t> nodes_;
  analyzers::utils::workload<func_t, max_num_steps> workload_;
  network::incomplete_signal_map<uint32_t, Ntk> signal_to_activity_;
  std::vector<analyzers::utils::signal_switching<func_t, max_num_steps>> activity_;
  analyzers::trackers::arrival_times_tracker<Ntk> arrival_;
  analyzers::trackers::sensing_times_tracker<Ntk> sensing_;
  analyzers::trackers::gate_load_tracker<Ntk> loading_;
  WinMngr& win_manager_;
};

} /* namespace profilers */

} /* namespace opto */

} /* namespace rinox */