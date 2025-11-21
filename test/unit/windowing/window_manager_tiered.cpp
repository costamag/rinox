#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <rinox/network/network.hpp>
#include <rinox/windowing/window_manager_tiered.hpp>

std::string const test_library =
    "GATE   inv1    1.0 O=!a ;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
    "GATE   and2    1.0 O=a*b;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
    "GATE   xor2    1.0 O=a^b;         PIN * INV 1   999 3.0 0.0 3.0 0.0";

struct WindowManagerParams1 : rinox::windowing::WindowManagerParamsDefault {
  static constexpr uint32_t kInpsCapacity = 32u;
  static constexpr uint32_t kDivsCapacity = 32u;
  static constexpr uint32_t kMffcCapacity = 32u;
  static constexpr uint32_t kTfosCapacity = 32u;
  static constexpr uint32_t kOutsCapacity = 32u;
  static constexpr uint32_t kMaxNumLeaves = 3u;
};

struct WindowManagerParams2 : rinox::windowing::WindowManagerParamsDefault {
  static constexpr uint32_t kInpsCapacity = 32u;
  static constexpr uint32_t kDivsCapacity = 32u;
  static constexpr uint32_t kMffcCapacity = 32u;
  static constexpr uint32_t kTfosCapacity = 32u;
  static constexpr uint32_t kOutsCapacity = 32u;
  static constexpr uint32_t kMaxNumLeaves = 8u;
};

TEST_CASE("Window made of only MFFC nodes and PIs - tiered version",
          "[window_manager_tiered]") {
  using Ntk =
      rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED,
                                    2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in(test_library);
  auto result = lorina::read_genlib(in, mockturtle::genlib_reader(gates));
  CHECK(result == lorina::return_code::success);

  Ntk ntk(gates);

  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  /* Database construction from network */
  fs.emplace_back();                                 // 0
  fs.emplace_back();                                 // 1
  fs.push_back(ntk.create_pi());                     // 2
  fs.push_back(ntk.create_pi());                     // 3
  fs.push_back(ntk.create_pi());                     // 4
  fs.push_back(ntk.create_node({fs[2], fs[3]}, 1));  // 5
  fs.push_back(ntk.create_node({fs[5], fs[4]}, 1));  // 6
  fs.push_back(ntk.create_node({fs[6], fs[5]}, 1));  // 7
  fs.push_back(ntk.create_node({fs[6], fs[7]}, 1));  // 8

  ntk.create_po(fs[8]);

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk(ntk);

  WindowManagerParams1 ps1;
  ps1.max_num_levels_tfo = 0u;
  rinox::windowing::WindowManagerStatsDefault st;

  rinox::windowing::WindowManager<DNtk, WindowManagerParams1> window(dntk, ps1,
                                                                     st);

  CHECK(window.Run(dntk.get_node(fs[8])));
  CHECK(window.IsValid());
  auto inps = window.GetSignalsInps();
  rinox::math::static_vector<typename Ntk::signal, 32u> inps_exp{8, 6, 4};
  CHECK(inps == inps_exp);
  auto mffc = window.GetNodesMffc();
  rinox::math::static_vector<typename Ntk::node, 32u> mffc_exp{5, 6, 7, 8};
  CHECK(mffc == mffc_exp);
}

TEST_CASE(
    "Window made of only MFFC nodes exceeding leaves capacity  - tiered "
    "version",
    "[window_manager_tiered]") {
  using Ntk =
      rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED,
                                    2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in(test_library);
  auto result = lorina::read_genlib(in, mockturtle::genlib_reader(gates));
  CHECK(result == lorina::return_code::success);

  Ntk ntk(gates);

  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  /* Database construction from network */
  fs.emplace_back();                                  // 0
  fs.emplace_back();                                  // 1
  fs.push_back(ntk.create_pi());                      // 2
  fs.push_back(ntk.create_pi());                      // 3
  fs.push_back(ntk.create_pi());                      // 4
  fs.push_back(ntk.create_pi());                      // 5
  fs.push_back(ntk.create_pi());                      // 6
  fs.push_back(ntk.create_node({fs[3]}, 0));          // 7
  fs.push_back(ntk.create_node({fs[2], fs[7]}, 1));   // 8
  fs.push_back(ntk.create_node({fs[4], fs[8]}, 1));   // 9
  fs.push_back(ntk.create_node({fs[8], fs[5]}, 1));   // 10
  fs.push_back(ntk.create_node({fs[6], fs[10]}, 1));  // 11
  fs.push_back(ntk.create_node({fs[9], fs[11]}, 1));  // 12

  ntk.create_po(fs[12]);

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk(ntk);

  WindowManagerParams1 ps1;
  ps1.max_num_levels_tfo = 0u;
  rinox::windowing::WindowManagerStatsDefault st;

  rinox::windowing::WindowManager<DNtk, WindowManagerParams1> window(dntk, ps1,
                                                                     st);

  CHECK(window.Run(dntk.get_node(fs[12])));
  CHECK(window.IsValid());
  auto const mffc = window.GetNodesMffc();
  CHECK(mffc == rinox::math::static_vector<typename Ntk::node, 32u>{9, 12});
  auto const inps = window.GetSignalsInps();
  CHECK(inps ==
        rinox::math::static_vector<typename Ntk::signal, 32u>{8, 16, 22});
}

TEST_CASE("Window where leaves expansion is possible - tiered version",
          "[window_manager_tiered]") {
  using Ntk =
      rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED,
                                    2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in(test_library);
  auto result = lorina::read_genlib(in, mockturtle::genlib_reader(gates));
  CHECK(result == lorina::return_code::success);

  Ntk ntk(gates);

  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  /* Database construction from network */
  fs.emplace_back();                                  // 0
  fs.emplace_back();                                  // 1
  fs.push_back(ntk.create_pi());                      // 2
  fs.push_back(ntk.create_pi());                      // 3
  fs.push_back(ntk.create_pi());                      // 4
  fs.push_back(ntk.create_pi());                      // 5
  fs.push_back(ntk.create_pi());                      // 6
  fs.push_back(ntk.create_node({fs[3]}, 0));          // 7
  fs.push_back(ntk.create_node({fs[2], fs[7]}, 1));   // 8
  fs.push_back(ntk.create_node({fs[4], fs[8]}, 1));   // 9
  fs.push_back(ntk.create_node({fs[8], fs[5]}, 1));   // 10
  fs.push_back(ntk.create_node({fs[6], fs[10]}, 1));  // 11
  fs.push_back(ntk.create_node({fs[9], fs[11]}, 1));  // 12

  ntk.create_po(fs[9]);
  ntk.create_po(fs[10]);
  ntk.create_po(fs[12]);

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::WindowManagerStatsDefault st;
  DNtk dntk(ntk);

  WindowManagerParams2 ps2;
  ps2.max_num_levels_tfo = 0u;

  rinox::windowing::WindowManager<DNtk, WindowManagerParams2> window(dntk, ps2,
                                                                     st);

  CHECK(window.Run(dntk.get_node(fs[12])));
  CHECK(window.IsValid());
  auto const mffc = window.GetNodesMffc();
  CHECK(mffc == rinox::math::static_vector<typename Ntk::node, 32u>{11, 12});
  auto const inps = window.GetSignalsInps();
  CHECK(inps ==
        rinox::math::static_vector<typename Ntk::signal, 32u>{6, 10, 8, 4, 12});
  auto const divs = window.GetSignalsDivs();
  CHECK(divs == rinox::math::static_vector<typename Ntk::signal, 32u>{
                    6, 10, 8, 4, 12, 14, 16, 20, 18});
}

TEST_CASE(
    "Window where leaves expansion is possible and side divs are present - "
    "tiered version",
    "[window_manager_tiered]") {
  using Ntk =
      rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED,
                                    2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in(test_library);
  auto result = lorina::read_genlib(in, mockturtle::genlib_reader(gates));
  CHECK(result == lorina::return_code::success);

  Ntk ntk(gates);

  using Signal = typename Ntk::signal;
  std::vector<Signal> fs;
  /* Database construction from network */
  fs.emplace_back();                                  // 0
  fs.emplace_back();                                  // 1
  fs.push_back(ntk.create_pi());                      // 2
  fs.push_back(ntk.create_pi());                      // 3
  fs.push_back(ntk.create_pi());                      // 4
  fs.push_back(ntk.create_pi());                      // 5
  fs.push_back(ntk.create_pi());                      // 6
  fs.push_back(ntk.create_node({fs[3]}, 0));          // 7
  fs.push_back(ntk.create_node({fs[2], fs[7]}, 1));   // 8
  fs.push_back(ntk.create_node({fs[4], fs[8]}, 1));   // 9
  fs.push_back(ntk.create_node({fs[8], fs[5]}, 1));   // 10
  fs.push_back(ntk.create_node({fs[6], fs[10]}, 1));  // 11
  fs.push_back(ntk.create_node({fs[9], fs[11]}, 1));  // 12
  fs.push_back(ntk.create_node({fs[9], fs[10]}, 1));  // 13
  fs.push_back(ntk.create_node({fs[13], fs[5]}, 1));  // 14

  ntk.create_po(fs[12]);
  ntk.create_po(fs[14]);

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::WindowManagerStatsDefault st;
  DNtk dntk(ntk);

  WindowManagerParams2 ps2;
  ps2.max_num_levels_tfo = 0u;

  rinox::windowing::WindowManager<DNtk, WindowManagerParams2> window(dntk, ps2,
                                                                     st);

  CHECK(window.Run(dntk.get_node(fs[12])));
  CHECK(window.IsValid());
  auto const mffc = window.GetNodesMffc();
  CHECK(mffc == rinox::math::static_vector<typename Ntk::node, 32u>{11, 12});
  auto const inps = window.GetSignalsInps();
  CHECK(inps ==
        rinox::math::static_vector<typename Ntk::signal, 32u>{6, 10, 8, 4, 12});
  auto const divs = window.GetSignalsDivs();
  CHECK(divs == rinox::math::static_vector<typename Ntk::signal, 32u>{
                    6, 10, 8, 4, 12, 14, 16, 20, 18, 26});
}
