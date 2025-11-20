#include <fmt/format.h>

#include <iostream>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
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
#include <rinox/opto/array/rewrite/array_rewrite_area.hpp>
#include <rinox/windowing/window_manager.hpp>
#include <string>
#include <vector>

#include "../../../../../common/baseline/mapping/lut_map.hpp"
#include "../../../../../common/baseline/preprocessing/aig_opto.hpp"
#include "../../../../../common/context.hpp"
#include "../../../../../common/env.hpp"
#include "../../../../../common/rinox/opto/algorithms/array/area_struct.hpp"

struct StaticParams : rinox::opto::array::AreaStructStaticParams {
  static constexpr uint32_t kMaxNumLeaves = 15u;
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr uint32_t kMaxInpsSolo = 6u;
};

rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 1u>
abc_lut_map(std::string const&,
            std::vector<rinox::experiments::aig_opto_t> const&,
            rinox::experiments::lut_map_t const&, std::string, int);

rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 1u>
abc_lutpack(std::string const&,
            std::vector<rinox::experiments::aig_opto_t> const&,
            rinox::experiments::lut_map_t const&, std::string, int);

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
  std::string const exp_path = std::string(RINOX_EXPERIMENTS_DIR) +
                               "/studies/resynthesis/array_based/area/struct/";
  std::string const config_path = exp_path + "config.json";

  int K0 = 6;
  std::string nameKLut = "|" + std::to_string(K0) + "LUT|abc";

  rinox::experiments::experiment<std::string, int, int, float, bool, int, int,
                                 float, bool>
      exp(exp_path, "report", "benchmark", "A(SOTA)", "D(SOTA)", "t(SOTA)",
          "E(SOTA)", "A(OPTO)", "D(OPTO)", "t(OPTO)", "E(OPTO)");

  auto doc = rinox::experiments::load_json_doc(config_path, &diag);
  if (!doc) return 1;

  auto ctx = rinox::experiments::load_context(*doc, &diag);
  auto const aig_opts = rinox::experiments::parse_aig_opto(*doc, &diag);
  auto map = rinox::experiments::load_lut_map_params(*doc, &diag);
  using DynamicParams = rinox::opto::array::AreaStructDynamicParams;
  auto rsy = rinox::experiments::array::load_area_struct_params<DynamicParams>(
      *doc, &diag);

  /* ==========================================================================
   */
  /*   Technology setup */
  /* ==========================================================================
   */
  auto static constexpr design_t = rinox::network::design_type_t::ARRAY_BASED;
  static constexpr uint32_t max_num_outputs = 1u;
  using Ntk = rinox::network::bound_network<design_t, max_num_outputs>;

  std::vector<std::array<double, max_num_outputs>> distr_bef;
  std::vector<std::array<double, max_num_outputs>> distr_aft;

  bool const success = rinox::experiments::foreach_benchmark(
      ctx, [&](std::string const& path, std::string const& name) {
        rinox::diagnostics::REPORT_DIAG(&diag, lorina::diagnostic_level::note,
                                        "Processing benchmark `{}`", name);

        /* ==========================================================================
         */
        /*   Baseline evaluation */
        /* ==========================================================================
         */
        std::cout << "mapping " << std::endl;
        Ntk ntk1 = abc_lutpack(path, aig_opts, map, name, K0);
        mockturtle::depth_view dntk1{ntk1};
        int depth1 = static_cast<int>(dntk1.depth());
        int area1 = static_cast<int>(ntk1.area());
        float time1 = 0.0;  // TODO: get time from abc_lutpack
        mockturtle::klut_network klut1 =
            rinox::network::convert_bound_to_klut<max_num_outputs>(ntk1);
        bool cec1 = rinox::experiments::abc_cec(klut1, path);
        std::cout << "SOTA: " << area1 << " " << depth1 << " " << time1 << " ";
        if (!cec1 && (name != "hyp")) std::cout << "NOT ";
        std::cout << "EQUIVALENT" << std::endl;
        /* ==========================================================================
         */
        /*   New algorithm evaluation */
        /* ==========================================================================
         */
        Ntk ntk2 = abc_lut_map(path, aig_opts, map, name, K0);
        rinox::opto::array::AreaStructStats rst;
        rinox::opto::array::AreaStruct<Ntk, StaticParams>(ntk2, rsy.ps, &rst);

        mockturtle::depth_view dntk2{ntk2};
        int depth2 = static_cast<int>(dntk2.depth());
        int area2 = static_cast<int>(ntk2.area());
        float time2 = 0.0;  // TODO: get time from abc_lutpack

        mockturtle::klut_network klut2 =
            rinox::network::convert_bound_to_klut<max_num_outputs>(ntk2);
        bool cec2 = rinox::experiments::abc_cec(klut2, path);
        std::cout << "OPTO: " << area2 << " " << depth2 << " " << time2 << " ";

        if (!cec2 && (name != "hyp")) std::cout << "NOT ";
        std::cout << "EQUIVALENT" << std::endl;
        /* ==========================================================================
         */
        /*   Results storage */
        /* ==========================================================================
         */
        exp(name, area1, depth1, time1, cec1, area2, depth2, time2, cec2);

        return true;
      });

  exp.save();
  exp.table();

  return 0;
}

rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 1u>
abc_lut_map(std::string const& path,
            std::vector<rinox::experiments::aig_opto_t> const& aig_optos,
            rinox::experiments::lut_map_t const& map, std::string name, int K) {
  std::string command = fmt::format(
      "abc -q \"read_aiger {}; strash; resyn; resyn2; if -K 6; write_blif "
      "/tmp/abc2_{}.blif\"",
      path, name);

#if WIN32
  std::unique_ptr<FILE, decltype(&_pclose)> pipe(_popen(command.c_str(), "r"),
                                                 _pclose);
#else
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(command.c_str(), "r"),
                                                pclose);
#endif

  if (!pipe) {
    throw std::runtime_error("popen() failed");
  }

  // Wait for ABC to finish (reading ensures process runs to completion)
  char buffer[128];
  while (fgets(buffer, sizeof(buffer), pipe.get())) {
    std::cout << buffer;  // optional: log ABC output
  }

  mockturtle::klut_network klut;
  std::string path_to_res = "/tmp/abc2_" + name + ".blif";
  if (lorina::read_blif(path_to_res, mockturtle::blif_reader(klut)) !=
      lorina::return_code::success) {
    throw std::runtime_error("Failed to read BLIF output from ABC");
  }
  std::cout << "#kLUT5 = " << klut.num_gates() << std::endl;
  return rinox::network::convert_klut_to_bound<
      rinox::network::design_type_t::ARRAY_BASED, 1u>(klut);
}

rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 1u>
abc_lutpack(std::string const& path,
            std::vector<rinox::experiments::aig_opto_t> const& aig_optos,
            rinox::experiments::lut_map_t const& map, std::string name, int K) {
  std::string command = fmt::format(
      "abc -q \"read_aiger {}; strash; resyn; resyn2; if -K 6; lutpack -S 3 -L "
      "100; lutpack -S 3 -L 100; write_blif "
      "/tmp/abc1_{}.blif\"",
      path, name);

#if WIN32
  std::unique_ptr<FILE, decltype(&_pclose)> pipe(_popen(command.c_str(), "r"),
                                                 _pclose);
#else
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(command.c_str(), "r"),
                                                pclose);
#endif

  if (!pipe) {
    throw std::runtime_error("popen() failed");
  }

  // Wait for ABC to finish (reading ensures process runs to completion)
  char buffer[128];
  while (fgets(buffer, sizeof(buffer), pipe.get())) {
    std::cout << buffer << std::endl;  // optional: log ABC output
  }

  mockturtle::klut_network klut;
  std::string path_to_res = "/tmp/abc1_" + name + ".blif";
  if (lorina::read_blif(path_to_res, mockturtle::blif_reader(klut)) !=
      lorina::return_code::success) {
    throw std::runtime_error("Failed to read BLIF output from ABC");
  }

  return rinox::network::convert_klut_to_bound<
      rinox::network::design_type_t::ARRAY_BASED, 1u>(klut);
}
