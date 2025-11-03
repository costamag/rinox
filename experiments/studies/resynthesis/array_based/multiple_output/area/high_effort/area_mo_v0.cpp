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
#include <rinox/opto/array/dual_output_resyn.hpp>
#include <rinox/windowing/window_manager.hpp>
#include <string>
#include <vector>

#include "../../../../../../common/baseline/mapping/lut_map.hpp"
#include "../../../../../../common/baseline/preprocessing/aig_opto.hpp"
#include "../../../../../../common/context.hpp"
#include "../../../../../../common/env.hpp"
#include "../../../../../../common/rinox/opto/algorithms/resynthesis.hpp"

template <uint32_t NumLeaves>
struct ResynParams
    : rinox::opto::algorithms::default_dual_output_resyn_params<NumLeaves> {
  static constexpr bool do_strashing = true;
  static constexpr uint32_t num_vars_sign = 12;
  static constexpr uint32_t max_cuts_size = 4u;
  static constexpr uint32_t max_cube_spfd = 10u;

  struct window_manager_params
      : rinox::windowing::default_window_manager_params {
    static constexpr uint32_t max_num_leaves = NumLeaves;
    bool preserve_depth = false;
    int32_t odc_levels = 0u;
    uint32_t skip_fanout_limit_for_divisors = 100u;
    uint32_t max_num_divisors{128};
  };

  window_manager_params window_manager_ps;
};


template<typename Ntk>
void RecordFaninsDistrib(Ntk&, std::vector<std::array<double, Ntk::kMaxNumOutputs>>&);

rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 2u>
abc_lut_map(std::string const&,
            std::vector<rinox::experiments::aig_opto_t> const&,
            rinox::experiments::lut_map_t const&, std::string, int);

struct perf_t {
  int nLuts;
  int depth;
  double fDual;
  double dArea;
  float rTime;
};

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
      "/studies/resynthesis/array_based/multiple_output/area/high_effort/";
  std::string const config_path = exp_path + "config.json";

  int K0 = 6;
  std::string nameKLut = "|" + std::to_string(K0) + "LUT|abc";

  rinox::experiments::experiment<std::string, int, int, int, int, double,
                                 double, float, int, int, double, double, float,
                                 bool, bool>
      exp(exp_path, "report", "benchmark", nameKLut, "|6LUT*|abc",
          "|65LUT|sota", "D(65LUT)|sota", "fDUAL|sota", "dA|sota", "t|sota",
          "|65LUT|opto", "D(65LUT)|opto", "fDUAL|opto", "dA|opto", "t|opto",
          "EQ(sota)", "EQ(opto)");

  auto doc = rinox::experiments::load_json_doc(config_path, &diag);
  if (!doc) return 1;

  auto ctx = rinox::experiments::load_context(*doc, &diag);
  auto const aig_opts = rinox::experiments::parse_aig_opto(*doc, &diag);
  auto map = rinox::experiments::load_lut_map_params(*doc, &diag);
  auto rsy =
      rinox::experiments::load_resynthesis_params<ResynParams<12>>(*doc, &diag);

  /* ==========================================================================
   */
  /*   Technology setup */
  /* ==========================================================================
   */
  auto static constexpr design_t = rinox::network::design_type_t::ARRAY_BASED;
  static constexpr uint32_t max_num_outputs = 2u;
  using Ntk = rinox::network::bound_network<design_t, max_num_outputs>;

  int nLutsAbc5 = 0, nLutsAbc6 = 0;
  int nLutsSota = 0, nLutsOpto = 0;
  int depthSota = 0, depthOpto = 0;
  double fDualSota = 0, fDualOpto = 0;
  double dAreaSota = 0, dAreaOpto = 0;
  float rTimeSota = 0, rTimeOpto = 0;
  double avgAreaOpto = 0;
  double avgAreaSota = 0;
  double avgDualOpto = 0;
  double avgDualSota = 0;
  double avgRunTimeOpto = 0;
  double avgRunTimeSota = 0;
  double avgDepthOpto = 0;
  double avgDepthSota = 0;
  double count = 0;

  std::vector<std::array<double, max_num_outputs>> distr_bef;
  std::vector<std::array<double, max_num_outputs>> distr_aft;

  std::unordered_map<std::string, perf_t> sota;
  sota["adder"] = {179, 51, 44.69, 30.35, 0.29};
  sota["arbiter"] = {2468, 18, 10.29, 9.33, 2.98};
  sota["bar"] = {448, 4, 14.29, 12.50, 0.39};
  sota["cavlc"] = {105, 4, 12.38, 11.02, 0.10};
  sota["ctrl"] = {18, 2, 61.11, 37.93, 0.06};
  sota["dec"] = {144, 2, 99.31, 49.83, 0.090};
  sota["div"] = {7279, 866, 38.27, 26.70, 56.31};
  sota["hyp"] = {24617, 4198, 81.58, 44.48, 399.65};
  sota["i2c"] = {242, 4, 30.58, 23.42, 0.21};
  sota["int2float"] = {41, 3, 14.63, 12.77, 0.07};
  sota["log2"] = {6476, 76, 27.53, 19.82, 16.26};
  sota["max"] = {620, 39, 28.55, 21.82, 0.54};
  sota["mem_ctrl"] = {9475, 30, 27.91, 19.34, 59.88};
  sota["multiplier"] = {4324, 54, 37.44, 26.98, 9.98};
  sota["priority"] = {159, 26, 12.58, 10.67, 0.14};
  sota["router"] = {72, 6, 5.56, 5.26, 0.07};
  sota["sin"] = {1203, 40, 25.35, 18.55, 1.33};
  sota["sqrt"] = {3908, 1118, 42.14, 15.87, 9.71};
  sota["square"] = {2692, 50, 46.84, 31.83, 6.34};
  sota["voter"] = {1380, 16, 28.48, 20.51, 0.72};

  std::unordered_map<std::string, int> abc_ifK6;
  abc_ifK6["adder"] = 257;
  abc_ifK6["arbiter"] = 2722;
  abc_ifK6["bar"] = 512;
  abc_ifK6["cavlc"] = 118;
  abc_ifK6["ctrl"] = 29;
  abc_ifK6["dec"] = 287;
  abc_ifK6["div"] = 9931;
  abc_ifK6["hyp"] = 0;
  abc_ifK6["i2c"] = 316;
  abc_ifK6["int2float"] = 47;
  abc_ifK6["log2"] = 8077;
  abc_ifK6["max"] = 793;
  abc_ifK6["mem_ctrl"] = 11470;
  abc_ifK6["multiplier"] = 5922;
  abc_ifK6["priority"] = 178;
  abc_ifK6["router"] = 50;
  abc_ifK6["sin"] = 1477;
  abc_ifK6["sqrt"] = 4645;
  abc_ifK6["square"] = 3949;
  abc_ifK6["voter"] = 1736;

  bool const success = rinox::experiments::foreach_benchmark(
      ctx, [&](std::string const& path, std::string const& name) {
        rinox::diagnostics::REPORT_DIAG(&diag, lorina::diagnostic_level::note,
                                        "Processing benchmark `{}`", name);

        /* ==========================================================================
         */
        /*   Baseline evaluation */
        /* ==========================================================================
         */
        Ntk ntk = abc_lut_map(path, aig_opts, map, name, K0);
        int nLutsAbc = static_cast<int>(ntk.area());
        std::set<std::string> no_cec = {"hyp", "div"};
        mockturtle::klut_network klut_sota = convert_bound_to_klut(ntk);
        bool const cec_sota = no_cec.find(name) == no_cec.end()
                                  ? rinox::experiments::abc_cec(klut_sota, path)
                                  : true;

        /* ==========================================================================
         */
        /*   New algorithm evaluation */
        /* ==========================================================================
         */
        RecordFaninsDistrib(ntk, distr_bef);

        rinox::opto::algorithms::resynthesis_stats rst;
        rinox::opto::algorithms::area_dual_output_resyn(ntk, rsy.ps, &rst);
        RecordFaninsDistrib(ntk, distr_aft);

        fDualOpto = 100 * rst.fDual;
        nLutsOpto = static_cast<int>(ntk.area());
        dAreaOpto = 100 * (1 - ((double)nLutsOpto) / ((double)abc_ifK6[name]));
        rTimeOpto = std::max(0.01, 0.001 * rst.time);
        mockturtle::depth_view dntk{ntk};
        depthOpto = static_cast<int>(dntk.depth());
        // rinox::analyzers::trackers::arrival_times_tracker times_track( dntk
        // );
        /* ==========================================================================
         */
        /*   Equivalence checking */
        /* ==========================================================================
         */

        //    ntk.foreach_gate( [&]( auto const& n ){
        //      std::cout << n << "=f(";
        //      ntk.foreach_fanin( n, [&]( auto const& fi, auto ii ){
        //        std::cout << fi.index << "." << fi.output << " ";
        //      } );
        //      std::cout << ") : ";
        //      auto tts = ntk.get_functions( n );
        //      for ( int i = 0; i < tts.size(); ++i )
        //      {
        //        std::cout << "  " << i << " ";
        //        kitty::print_binary( tts[i] );
        //        std::cout << std::endl;
        //      }
        //    } );

        mockturtle::klut_network klut_opto =
            rinox::network::convert_bound_to_klut<max_num_outputs>(ntk);
        bool const cec_opto = no_cec.find(name) == no_cec.end()
                                  ? rinox::experiments::abc_cec(klut_opto, path)
                                  : true;

        //    std::cout << std::endl;
        //    std::cout << std::endl;
        //
        //    klut_opto.foreach_gate( [&]( auto const& n ){
        //      std::cout << n << "=f(";
        //      klut_opto.foreach_fanin( n, [&]( auto const& fi, auto ii ){
        //        std::cout << fi << " ";
        //      } );
        //      std::cout << ") : ";
        //      auto tt = klut_opto.node_function( n );
        //        std::cout << "  0 ";
        //        kitty::print_binary( tt );
        //        std::cout << std::endl;
        //      } );

        nLutsAbc6 = abc_ifK6[name];
        nLutsSota = sota[name].nLuts;
        depthSota = sota[name].depth;
        fDualSota = sota[name].fDual;
        dAreaSota = sota[name].dArea;
        rTimeSota = sota[name].rTime;
        /* ==========================================================================
         */
        /*   Results storage */
        /* ==========================================================================
         */
        exp(name, nLutsAbc, abc_ifK6[name], nLutsSota, depthSota, fDualSota,
            dAreaSota, rTimeSota, nLutsOpto, depthOpto, fDualOpto, dAreaOpto,
            rTimeOpto, cec_sota, cec_opto);

        avgAreaOpto += dAreaOpto;
        avgAreaSota += dAreaSota;
        avgDualOpto += fDualOpto;
        avgDualSota += fDualSota;
        avgRunTimeSota += rTimeSota;
        avgRunTimeOpto += rTimeOpto;
        count += 1;

        return true;
      });

  exp.save();
  exp.table();

  avgAreaOpto /= count;
  avgAreaSota /= count;
  avgDualOpto /= count;
  avgDualSota /= count;
  avgRunTimeSota /= count;
  avgRunTimeOpto /= count;

  std::cout << "\n==== Average Results over " << count << " benchmarks ====\n";

  std::cout << std::left << std::setw(20) << "Metric" << std::setw(15) << "SOTA"
            << std::setw(15) << "Proposed" << std::setw(15) << "Improvement"
            << "\n";

  std::cout << std::string(65, '-') << "\n";

  auto pct = [](double base, double newval) {
    return (base - newval) / base * 100.0;
  };

  std::cout << std::setw(20) << "Area" << std::setw(15) << avgAreaSota
            << std::setw(15) << avgAreaOpto << std::setw(14)
            << -pct(avgAreaSota, avgAreaOpto) << "%\n";

  std::cout << std::setw(20) << "Dual" << std::setw(15) << avgDualSota
            << std::setw(15) << avgDualOpto << std::setw(14)
            << -pct(avgDualSota, avgDualOpto) << "%\n";

  std::cout << std::setw(20) << "Runtime (s)" << std::setw(15) << avgRunTimeSota
            << std::setw(15) << avgRunTimeOpto << std::setw(14)
            << pct(avgRunTimeSota, avgRunTimeOpto) << "%\n";

  return 0;
}

rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 2u>
abc_lut_map(std::string const& path,
            std::vector<rinox::experiments::aig_opto_t> const& aig_optos,
            rinox::experiments::lut_map_t const& map, std::string name, int K) {
  std::string command = fmt::format(
      "abc -q \"read_aiger {}; resyn; resyn2; if -a -K {}; write_blif "
      "/tmp/abc_{}.blif\"",
      path, K, name);

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
  std::string path_to_res = "/tmp/abc_" + name + ".blif";
  if (lorina::read_blif(path_to_res, mockturtle::blif_reader(klut)) !=
      lorina::return_code::success) {
    throw std::runtime_error("Failed to read BLIF output from ABC");
  }

  std::cout << "#kLUT5 = " << klut.num_gates() << std::endl;
  using Ntk =
      rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED,
                                    2u>;

  return rinox::network::convert_klut_to_bound<
      rinox::network::design_type_t::ARRAY_BASED, 2u>(klut);
}

 template<typename Ntk>
  void RecordFaninsDistrib(Ntk&ntk, std::vector<std::array<double, Ntk::kMaxNumOutputs>>& distr)
  {
    std::array<double, Ntk::kMaxNumOutputs> zero;
    for(uint32_t i=0; i<Ntk::kMaxNumOutputs; ++i) zero[i]=0;

    double total = 0;
    ntk.foreach_gate([&](auto const& n){
      if(ntk.is_dead(n)) return true;
      uint32_t num_fanins = ntk.fanin_size(n);
      while (distr.size()<=num_fanins)
        distr.push_back(zero);
      distr[num_fanins][ntk.num_outputs(n)-1]+=1;
      total+=1;
      return true;
    });
    for(int num_fanins=1; num_fanins<distr.size();++num_fanins)
    {
      std::cout << "|FIs|: " << num_fanins << " ";
      for(int i=0; i<Ntk::kMaxNumOutputs;++i)
        std::cout << distr[num_fanins][i] << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }