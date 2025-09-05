#include "../common/baseline/mapping/tech_map.hpp"
#include "../common/baseline/preprocessing/aig_opto.hpp"
#include "../common/context.hpp"
#include "../common/env.hpp"
#include "../common/rinox/opto/algorithms/resynthesis.hpp"
#include <iostream>
#include <rinox/analyzers/trackers/arrival_times_tracker.hpp>
#include <rinox/databases/mapped_database.hpp>
#include <rinox/network/converters.hpp>
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

template<uint32_t NumLeaves>
struct ResynParams : rinox::opto::algorithms::default_resynthesis_params<NumLeaves>
{
  static constexpr bool do_strashing = true;
  static constexpr uint32_t num_vars_sign = 12;
  static constexpr uint32_t max_cuts_size = 4u;
  static constexpr uint32_t max_cube_spfd = 10u;

  struct window_manager_params : rinox::windowing::default_window_manager_params
  {
    static constexpr uint32_t max_num_leaves = NumLeaves;
    bool preserve_depth = false;
    int32_t odc_levels = 0u;
    uint32_t skip_fanout_limit_for_divisors = 100u;
    uint32_t max_num_divisors{ 128 };
  };

  window_manager_params window_manager_ps;
};

template<typename TechLib>
std::optional<std::tuple<mockturtle::cell_view<mockturtle::block_network>, mockturtle::emap_stats>> evaluate_sota( std::string const&, std::vector<rinox::experiments::aig_opto_t> const&, rinox::experiments::tech_map_t const&, TechLib const& );

int main()
{
  /* ========================================================================== */
  /*   Diagnostics setup                                                        */
  /* ========================================================================== */
  rinox::diagnostics::text_diagnostics client;
  lorina::diagnostic_engine diag( &client );

  /* ========================================================================== */
  /*   Experimental setup                                                       */
  /* ========================================================================== */
  std::string const exp_path = std::string( RINOX_EXPERIMENTS_DIR ) + "/studies/resynthesis/area/high_effort/";
  std::string const config_path = exp_path + "config.json";

  rinox::experiments::experiment<
      std::string,
      double, double,
      double, double,
      float, float,
      bool, bool>
      exp(
          exp_path,
          "report",
          "benchmark",
          "Area(sota)", "Area(opto)",
          "Delay(sota)", "Delay(opto)",
          "Time(sota)", "Time(opto)",
          "EQ(sota)", "EQ(opto)" );

  auto doc = rinox::experiments::load_json_doc( config_path, &diag );
  if ( !doc )
    return 1;

  auto ctx = rinox::experiments::load_context( *doc, &diag );
  auto const aig_opts = rinox::experiments::parse_aig_opto( *doc, &diag );
  auto map = rinox::experiments::load_tech_map_params( *doc, &diag );
  auto rsy = rinox::experiments::load_resynthesis_params<ResynParams<12>>( *doc, &diag );

  /* ========================================================================== */
  /*   Technology setup                                                         */
  /* ========================================================================== */
  rinox::diagnostics::REPORT_DIAG_RAW( &diag, lorina::diagnostic_level::note,
                                       "Processing technology library" );
  auto const& gates = rinox::experiments::load_gates( ctx, &diag );
  if ( !gates )
    return false;
  mockturtle::tech_library<9> tech_lib( *gates, map.tps );

  auto static constexpr design_t = rinox::network::design_type_t::CELL_BASED;
  static constexpr uint32_t max_num_outputs = 2u;
  using Ntk = rinox::network::bound_network<design_t, max_num_outputs>;
  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( *gates );
  rinox::databases::mapped_database<Ntk, 4u> db( lib );
  db.load( rsy.database.path );

  bool const success = rinox::experiments::foreach_benchmark( ctx, [&]( std::string const& path, std::string const& name ) {
    rinox::diagnostics::REPORT_DIAG( &diag, lorina::diagnostic_level::note,
                                     "Processing benchmark `{}`", name );

    /* ========================================================================== */
    /*   Baseline evaluation                                                      */
    /* ========================================================================== */
    auto sota = evaluate_sota( path, aig_opts, map, tech_lib );
    if ( !sota )
    return false;
    auto [ntk_sota, st_sota] = *sota;
    mockturtle::write_verilog_with_cell( ntk_sota, name + ".sota.v" );

    /* ========================================================================== */
    /*   New algorithm evaluation                                                 */
    /* ========================================================================== */
    Ntk ntk_opto = rinox::network::convert_mapped_to_bound( ntk_sota, *gates );
    mockturtle::depth_view dntk_opto{ ntk_opto };
    rinox::analyzers::trackers::arrival_times_tracker times_track( dntk_opto );

    rinox::opto::algorithms::resynthesis_stats rst;

    rinox::opto::algorithms::area_resynthesize( dntk_opto, db, rsy.ps, &rst );
    rinox::io::verilog::write_verilog( ntk_opto, name + ".opto.v" );

    /* ========================================================================== */
    /*   Equivalence checking                                                     */
    /* ========================================================================== */
    bool const cec_sota = rinox::experiments::abc_cec_mapped_cell( ntk_sota, path, ctx.technology.path );
    bool const cec_opto = rinox::experiments::abc_cec_rinox( ntk_opto, path, ctx.technology.path );

    /* ========================================================================== */
    /*   Results storage                                                          */
    /* ========================================================================== */
    exp( name,
         ntk_sota.compute_area(), ntk_opto.area(),
         ntk_sota.compute_worst_delay(), times_track.worst_delay(),
         mockturtle::to_seconds( st_sota.time_total ), 0,
         cec_sota, cec_opto );
    return true;
  } );

  exp.save();
  exp.table();
  return 0;
}

template<typename TechLib>
std::optional<std::tuple<mockturtle::cell_view<mockturtle::block_network>, mockturtle::emap_stats>> evaluate_sota( std::string const& path, std::vector<rinox::experiments::aig_opto_t> const& aig_optos, rinox::experiments::tech_map_t const& map, TechLib const& tech_lib )
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
