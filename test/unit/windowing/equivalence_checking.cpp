#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <rinox/network/network.hpp>
#include <rinox/windowing/equivalence_checking.hpp>
#include <rinox/windowing/window_manager.hpp>
#include <rinox/windowing/window_simulator.hpp>

std::string const test_library =
    "GATE   inv1    1.0 O=!a ;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
    "GATE   and2    1.0 O=a*b;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
    "GATE   or2     1.0 O=a+b;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
    "GATE   xor2    1.0 O=a^b;         PIN * INV 1   999 3.0 0.0 3.0 0.0\n"
    "GATE   nor2    1.0 O=!(a+b);      PIN * INV 1   999 3.0 0.0 3.0 0.0";

struct window_manager_params : rinox::windowing::default_window_manager_params {
  static constexpr uint32_t max_num_leaves = 8u;
};

TEST_CASE("CEC a small window", "[window_checking]") {
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
  fs.emplace_back(ntk.create_pi());  // 0
  fs.emplace_back(ntk.create_pi());  // 1

  fs.emplace_back(ntk.create_node({fs[0]}, 0));         // 2
  fs.emplace_back(ntk.create_node({fs[1]}, 0));         // 3
  fs.emplace_back(ntk.create_node({fs[2], fs[3]}, 1));  // 4
  fs.emplace_back(ntk.create_node({fs[4]}, 0));         // 5

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk(ntk);

  window_manager_params ps;

  using Manager = rinox::windowing::window_manager<DNtk>;
  Manager window(dntk, ps, st);
  CHECK(window.run(dntk.get_node(fs[5])));
  rinox::windowing::window_simulator sim(dntk);
  sim.run(window);

  rinox::windowing::CECResynParams cec_ps;
  rinox::windowing::CECResynStats cec_st;
  rinox::windowing::CECResyn<Ntk, Manager> cec(ntk, window, cec_ps, cec_st);

  fs.emplace_back(ntk.create_node({fs[0], fs[1]}, 2));  // 6
  std::vector<signal> other{fs[6]};
  auto res = cec.SolveSignals(other);
  CHECK(res);
  CHECK(*res);

  fs.emplace_back(ntk.create_node({fs[0], fs[1]}, 3));  // 7
  other = {fs[7]};
  auto res7 = cec.SolveSignals(other);
  CHECK(res7);
  CHECK(!(*res7));

  fs.emplace_back(ntk.create_node({fs[0], fs[1]}, 4));  // 8
  fs.emplace_back(ntk.create_node({fs[8]}, 0));         // 9
  other = {fs[9]};
  auto res9 = cec.SolveSignals(other);
  CHECK(res9);
  CHECK(*res9);
}

TEST_CASE("CEC a small window with ODCs", "[window_checking]") {
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
  fs.emplace_back(ntk.create_pi());  // 0
  fs.emplace_back(ntk.create_pi());  // 1

  fs.emplace_back(ntk.create_node({fs[0], fs[1]}, 1));  // 2
  fs.emplace_back(ntk.create_node({fs[0], fs[1]}, 3));  // 3
  fs.emplace_back(ntk.create_node({fs[2], fs[3]}, 2));  // 4
  ntk.create_po(fs[4]);

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk(ntk);

  window_manager_params ps;
  ps.odc_levels = 2;

  using Manager = rinox::windowing::window_manager<DNtk>;
  Manager window(dntk, ps, st);
  CHECK(window.run(dntk.get_node(fs[3])));
  rinox::windowing::window_simulator sim(dntk);
  sim.run(window);

  rinox::windowing::CECResynParams cec_ps;
  rinox::windowing::CECResynStats cec_st;
  rinox::windowing::CECResyn<Ntk, Manager> cec(ntk, window, cec_ps, cec_st);

  fs.emplace_back(ntk.create_node({fs[0], fs[1]}, 2));  // 5
  std::vector<signal> other{fs[5]};
  auto res = cec.SolveSignals(other);
  CHECK(res);
  CHECK(*res);
}

TEST_CASE("CEC a small two-outputs window with ODCs", "[window_checking]") {
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
  fs.emplace_back(ntk.create_pi());  // 0
  fs.emplace_back(ntk.create_pi());  // 1

  fs.emplace_back(ntk.create_node({fs[0], fs[1]}, 2));  // 2
  fs.emplace_back(ntk.create_node({fs[0], fs[1]}, 1));  // 3
  fs.emplace_back(ntk.create_node({fs[3]}, 0));  // 4
  fs.emplace_back(ntk.create_node({fs[2]}, 0));  // 5
  fs.emplace_back(ntk.create_node({fs[3], fs[5]}, 2));  // 6
  fs.emplace_back(ntk.create_node({fs[2], fs[4]}, 1));  // 7

  ntk.create_po(fs[6]);
  ntk.create_po(fs[7]);

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk(ntk);

  window_manager_params ps;
  ps.odc_levels = 2;

  using Manager = rinox::windowing::window_manager<DNtk>;
  Manager window(dntk, ps, st);
  CHECK(window.run(dntk.get_node(fs[2])));
  rinox::windowing::window_simulator sim(dntk);
  sim.run(window);

  rinox::windowing::CECResynParams cec_ps;
  rinox::windowing::CECResynStats cec_st;
  rinox::windowing::CECResyn<Ntk, Manager> cec(ntk, window, cec_ps, cec_st);

  fs.emplace_back(ntk.create_node({fs[0], fs[1]}, 3));  // 8
  std::vector<signal> other{fs[8]};
  auto res = cec.SolveSignals(other);
  CHECK(res);
  CHECK(*res);
}