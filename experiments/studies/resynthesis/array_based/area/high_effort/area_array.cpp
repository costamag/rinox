#include <fmt/format.h>

#include <iostream>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/utils/name_utils.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/names_view.hpp>
#include <rinox/analyzers/trackers/arrival_times_tracker.hpp>
#include <rinox/databases/mapped_database.hpp>
#include <rinox/network/converters.hpp>
#include <rinox/opto/algorithms/array/area_resyn.hpp>
#include <rinox/opto/algorithms/resynthesize.hpp>
#include <rinox/synthesis/lut/synth.hpp>
#include <rinox/windowing/window_manager.hpp>
#include <string>
#include <vector>

#include "../../../../../common/baseline/mapping/lut_map.hpp"
#include "../../../../../common/baseline/preprocessing/aig_opto.hpp"
#include "../../../../../common/context.hpp"
#include "../../../../../common/env.hpp"
#include "../../../../../common/rinox/opto/algorithms/array/area_resyn.hpp"

struct ResynStaticParams
    : rinox::opto::algorithms::array::DefaultAreaResynStaticParams {
  // static parameters
  static constexpr uint32_t kMaxNumOutputs = 1u;
  static constexpr uint32_t kMaxNumLeaves = 12u;
  static constexpr uint32_t kNumVarsSign = 12u;
  static constexpr uint32_t kMaxCutsSize = 8u;
  static constexpr uint32_t kMaxCubeSpfd = 6u;
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr uint32_t kMaxInpsSolo = 6u;
  static constexpr rinox::synthesis::lut::LutSynthMethod kLutSynth =
      rinox::synthesis::lut::LutSynthMethod::kShannonSolo;
};

std::optional<std::tuple<mockturtle::klut_network, mockturtle::lut_map_stats>>
evaluate_sota(std::string const&,
              std::vector<rinox::experiments::aig_opto_t> const&,
              rinox::experiments::lut_map_t const&);

int main() {
  /* ==========================================================================
   */
  /*   Diagnostics setup */
  /* ==========================================================================
   */
  rinox::diagnostics::text_diagnostics client;
  lorina::diagnostic_engine diag(&client);

  /* ==========================================================================
   */
  /*   Experimental setup */
  /* ==========================================================================
   */
  std::string const exp_path =
      std::string(RINOX_EXPERIMENTS_DIR) +
      "/studies/resynthesis/array_based/area/high_effort/";
  std::string const config_path = exp_path + "config.json";

  rinox::experiments::experiment<std::string, double, double, double, double,
                                 float, float, bool, bool>
      exp(exp_path, "report", "benchmark", "Area(sota)", "Area(opto)",
          "Delay(sota)", "Delay(opto)", "Time(sota)", "Time(opto)", "EQ(sota)",
          "EQ(opto)");

  auto doc = rinox::experiments::load_json_doc(config_path, &diag);
  if (!doc) return 1;

  auto ctx = rinox::experiments::load_context(*doc, &diag);
  auto const aig_opts = rinox::experiments::parse_aig_opto(*doc, &diag);
  auto map = rinox::experiments::load_lut_map_params(*doc, &diag);
  using AreaResynParams = rinox::opto::algorithms::array::AreaResynParams;
  auto rsy =
      rinox::experiments::array::load_area_resyn_params<AreaResynParams>(
          *doc, &diag);

  /* ==========================================================================
   */
  /*   Technology setup */
  /* ==========================================================================
   */
  auto static constexpr design_t = rinox::network::design_type_t::ARRAY_BASED;
  using Ntk = rinox::network::bound_network<design_t,
                                            ResynStaticParams::kMaxNumOutputs>;

  bool const success = rinox::experiments::foreach_benchmark(
      ctx, [&](std::string const& path, std::string const& name) {
        rinox::diagnostics::REPORT_DIAG(&diag, lorina::diagnostic_level::note,
                                        "Processing benchmark `{}`", name);

        /* ==========================================================================
         */
        /*   Baseline evaluation */
        /* ==========================================================================
         */
        auto sota = evaluate_sota(path, aig_opts, map);
        if (!sota) return false;
        auto [klut_sota, st_sota] = *sota;

        /* ==========================================================================
         */
        /*   New algorithm evaluation */
        /* ==========================================================================
         */
        Ntk ntk_opto = rinox::network::convert_klut_to_bound<
            design_t, ResynStaticParams::kMaxNumOutputs>(klut_sota);
        double const area_before = ntk_opto.area();
        rinox::opto::algorithms::array::AreaResynStats rst;
        rinox::opto::algorithms::array::AreaResyn<Ntk, ResynStaticParams>(
            ntk_opto, rsy.ps, &rst);
        mockturtle::depth_view dntk_opto{ntk_opto};
        rinox::analyzers::trackers::arrival_times_tracker times_track(
            dntk_opto);
        rst.Report();
        /* ==========================================================================
         */
        /*   Equivalence checking */
        /* ==========================================================================
         */
        mockturtle::klut_network klut_opto =
            rinox::network::convert_bound_to_klut<
                ResynStaticParams::kMaxNumOutputs>(ntk_opto);
        bool const cec_sota = rinox::experiments::abc_cec(klut_sota, path);
        bool const cec_opto = rinox::experiments::abc_cec(klut_opto, path);

        /* ==========================================================================
         */
        /*   Results storage */
        /* ==========================================================================
         */
        exp(name, area_before, ntk_opto.area(), st_sota.delay,
            times_track.worst_delay(),
            mockturtle::to_seconds(st_sota.time_total), rst.GetTimeTotal(),
            cec_sota, cec_opto);
        return true;
      });

  exp.save();
  exp.table();
  return 0;
}

std::optional<std::tuple<mockturtle::klut_network, mockturtle::lut_map_stats>>
evaluate_sota(std::string const& path,
              std::vector<rinox::experiments::aig_opto_t> const& aig_optos,
              rinox::experiments::lut_map_t const& map) {
  // load the AIG
  mockturtle::names_view<mockturtle::aig_network> aig;
  if (lorina::read_aiger(path, mockturtle::aiger_reader(aig)) !=
      lorina::return_code::success) {
    return std::nullopt;
  }

  // preprocess the AIG
  for (auto const& aig_opt : aig_optos)
    rinox::experiments::apply_aig_opto(aig, aig_opt);

  // map to technology
  mockturtle::lut_map_stats mst;
  mockturtle::klut_network mapped = mockturtle::lut_map(aig, map.mps, &mst);
  return std::make_tuple(mapped, mst);
}
