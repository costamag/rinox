#include "../common/baseline/mapping/tech_map.hpp"
#include "../common/baseline/preprocessing/aig_opto.hpp"
#include "../common/context.hpp"
#include "../common/env.hpp"
#include "../common/rinox/windowing/windowing.hpp"
#include <iostream>
#include <rinox/analyzers/trackers/arrival_times_tracker.hpp>
#include <rinox/network/converters.hpp>
#include <rinox/databases/mapped_database.hpp>
#include <rinox/opto/algorithms/resynthesize.hpp>
#include <rinox/windowing/window_manager.hpp>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/utils/name_utils.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/cell_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/names_view.hpp>
#include <mockturtle/algorithms/resubstitution.hpp>

/* ========================================================================== */
/*   Define the stats for the experiment                                      */
/* ========================================================================== */
struct windowing_stats_t
{
  uint32_t count = 0;
  double avg_num_leaves = 0;
  double avg_num_divs = 0;
  double avg_num_mffcs = 0;
  float runtime = 0;
};

template<uint32_t NumLeaves>
struct WindowParams : rinox::windowing::default_window_manager_params
{
  static constexpr uint32_t max_num_leaves = NumLeaves;
  bool preserve_depth = false;
  int32_t odc_levels = 0u;
  uint32_t skip_fanout_limit_for_divisors = 100u;
  uint32_t max_num_divisors{ 128 };
};

template<typename TechLib>
std::optional<std::tuple<mockturtle::cell_view<mockturtle::block_network>, mockturtle::emap_stats>> init_network( std::string const&, std::vector<rinox::experiments::aig_opto_t> const&, rinox::experiments::tech_map_t const&, TechLib const& );

template<typename Ntk>
void sota_analysis( Ntk &, mockturtle::resubstitution_params const&, windowing_stats_t & );

template<typename Ntk, typename Params>
void opto_analysis( Ntk &, Params const&, windowing_stats_t & );

int main()
{
  static constexpr uint32_t MaxNumLeaves = 12u;
  /* ========================================================================== */
  /*   Diagnostics setup                                                        */
  /* ========================================================================== */
  rinox::diagnostics::text_diagnostics client;
  lorina::diagnostic_engine diag( &client );

  /* ========================================================================== */
  /*   Experimental setup                                                       */
  /* ========================================================================== */
  std::string const exp_path = std::string( RINOX_EXPERIMENTS_DIR ) + "/studies/windowing/";
  std::string const config_path = exp_path + "config.json";

  rinox::experiments::experiment<
      std::string,
      uint32_t, uint32_t,
      double, double,
      double, double,
      double, double,
      float, float>
      exp(
          exp_path,
          "windowing_report",
          "benchmark",
          "N(sota)", "N(opto)",
          "<L(sota)>", "<L(opto)>",
          "<D(sota)>", "<D(opto)>",
          "<M(sota)>", "<M(opto)>",
          "T(sota)", "T(opto)" );

  auto doc = rinox::experiments::load_json_doc( config_path, &diag );
  if ( !doc )
    return 1;

  auto ctx = rinox::experiments::load_context( *doc, &diag );
  auto const aig_opts = rinox::experiments::parse_aig_opto( *doc, &diag );
  auto map = rinox::experiments::load_tech_map_params( *doc, &diag );
  auto wps = rinox::experiments::load_windowing_params<MaxNumLeaves>( *doc, &diag );

  /* ========================================================================== */
  /*   Technology setup                                                         */
  /* ========================================================================== */
  rinox::diagnostics::REPORT_DIAG_RAW( &diag, lorina::diagnostic_level::note,
    "Processing technology library" );
  auto const& gates = rinox::experiments::load_gates( ctx, &diag );
  if ( !gates )
    return false;
  mockturtle::tech_library<9> tech_lib( *gates, map.tps );
  
  /* ========================================================================== */
  /*   SOTA parameters                                                          */
  /* ========================================================================== */
  mockturtle::resubstitution_params ps_sota;
  ps_sota.max_pis = MaxNumLeaves;
  ps_sota.preserve_depth = wps.preserve_depth;
  ps_sota.odc_levels = wps.odc_levels;
  ps_sota.skip_fanout_limit_for_divisors = wps.skip_fanout_limit_for_divisors;
  ps_sota.max_divisors = wps.max_num_divisors;
  
  /* ========================================================================== */
  /*   OPTO parameters                                                          */
  /* ========================================================================== */
  WindowParams<12u> ps_opto;
  ps_opto.preserve_depth = wps.preserve_depth;
  ps_opto.odc_levels = wps.odc_levels;
  ps_opto.skip_fanout_limit_for_divisors = wps.skip_fanout_limit_for_divisors;
  ps_opto.max_num_divisors = wps.max_num_divisors;
  ps_opto.max_reconvergence = wps.max_reconvergence;
  
  /* ========================================================================== */
  /*   Iterate over the benchmarks                                              */
  /* ========================================================================== */
  
  bool const success = rinox::experiments::foreach_benchmark( ctx, [&]( std::string const& path, std::string const& name ) {

  rinox::diagnostics::REPORT_DIAG( &diag, lorina::diagnostic_level::note,
    "Processing benchmark `{}`", name );

    /* ========================================================================== */
    /*   Network construction                                                     */
    /* ========================================================================== */
    auto init = init_network( path, aig_opts, map, tech_lib );
    if ( !init ) return false;
    auto [ntk_init, st_init] = *init;
    auto static constexpr design_t = rinox::network::design_type_t::CELL_BASED;
    static constexpr uint32_t max_num_outputs = 2u;
    using Ntk = rinox::network::bound_network<design_t, max_num_outputs>;
    Ntk ntk_opto = rinox::network::convert_mapped_to_bound<mockturtle::cell_view<mockturtle::block_network>, design_t, max_num_outputs>( ntk_init, *gates );
    mockturtle::depth_view dntk {ntk_opto};

    /* ========================================================================== */
    /*   SOTA analysis                                                            */
    /* ========================================================================== */
    windowing_stats_t st_sota;
    sota_analysis( dntk, ps_sota, st_sota );
    
    /* ========================================================================== */
    /*   OPTO analysis                                                            */
    /* ========================================================================== */
    windowing_stats_t st_opto;
    opto_analysis( dntk, ps_opto, st_opto );

    /* ========================================================================== */
    /*   Results storage                                                          */
    /* ========================================================================== */
    exp( name,
      st_sota.count, st_opto.count,
      st_sota.avg_num_leaves, st_opto.avg_num_leaves,
      st_sota.avg_num_divs, st_opto.avg_num_divs,
      st_sota.avg_num_mffcs, st_opto.avg_num_mffcs,
      st_sota.runtime, st_opto.runtime );
    return true;
  } );

  exp.save();
  exp.table();

  std::cout << "================================================" << std::endl;
  std::cout << "=  Legend:                                     =" << std::endl;
  std::cout << "================================================" << std::endl;
  std::cout << "=    N      : number of windows                =" << std::endl;
  std::cout << "=    T      : runtion                          =" << std::endl;
  std::cout << "=    <L>    : average number of leaves         =" << std::endl;
  std::cout << "=    <D>    : average number of divisors       =" << std::endl;
  std::cout << "=    <M>    : average number of MFFC nodes     =" << std::endl;
  std::cout << "=    (sota) : state-of-the-art                 =" << std::endl;
  std::cout << "=    (opto) : optimization under investigation =" << std::endl;
  std::cout << "================================================" << std::endl;

  return 0;
}

template<typename TechLib>
std::optional<std::tuple<mockturtle::cell_view<mockturtle::block_network>, mockturtle::emap_stats>> init_network( std::string const& path, std::vector<rinox::experiments::aig_opto_t> const& aig_optos, rinox::experiments::tech_map_t const& map, TechLib const& tech_lib )
{
  // load the AIG
  mockturtle::names_view<mockturtle::aig_network> aig;
  if ( lorina::read_aiger( path, mockturtle::aiger_reader( aig ) ) != lorina::return_code::success )
  {
    return std::nullopt;
  }

  // preprocess the AIG
  for ( auto const& aig_opt : aig_optos )
    rinox::experiments::apply_aig_opto( aig, aig_opt );

  // map to technology
  mockturtle::emap_stats mst;
  mockturtle::cell_view<mockturtle::block_network> mapped = mockturtle::emap<9>( aig, tech_lib, map.mps, &mst );
  return std::make_tuple( mapped, mst );
}


template<typename Ntk>
void sota_analysis( Ntk & ntk, mockturtle::resubstitution_params const& ps, windowing_stats_t & st )
{
  mockturtle::detail::default_collector_stats est;
  mockturtle::detail::default_divisor_collector win_mngr( ntk, ps, est );
  uint32_t mffc_res;
  ntk.foreach_gate( [&]( auto const& n ){
    if ( win_mngr.run( n, mffc_res ) )
    {
      auto const N = static_cast<double>(++st.count);

      st.avg_num_leaves = ( st.avg_num_leaves * ( N - 1 ) + win_mngr.leaves.size() ) / N;
      st.avg_num_divs = ( st.avg_num_divs * ( N - 1 ) + win_mngr.divs.size() ) / N;
      st.avg_num_mffcs = ( st.avg_num_mffcs * ( N - 1 ) + win_mngr.mffc.size() ) / N;
    }
    
    return true;
  } );
}

template<typename Ntk, typename Params>
void opto_analysis( Ntk & ntk, Params const& ps, windowing_stats_t & st )
{
  rinox::windowing::window_manager_stats est;
  rinox::windowing::window_manager win_mngr( ntk, ps, est );
  ntk.foreach_gate( [&]( auto const& n ){
    if ( win_mngr.run( n ) )
    {
      auto const N = static_cast<double>(++st.count);
      st.avg_num_leaves = ( st.avg_num_leaves * ( N - 1 ) + win_mngr.num_inputs() ) / N;
      st.avg_num_divs = ( st.avg_num_divs * ( N - 1 ) + win_mngr.num_divisors() ) / N;
      st.avg_num_mffcs = ( st.avg_num_mffcs * ( N - 1 ) + win_mngr.num_mffc() ) / N;
    }  
    return true;
  } );
}