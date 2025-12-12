#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <rinox/dependency/rewire_dependencies.hpp>
#include <rinox/network/network.hpp>
#include <rinox/windowing/window_manager.hpp>
#include <rinox/windowing/window_simulator.hpp>

std::string const test_library =
    "GATE   and2    1.0 O=a*b;                 PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   or2     1.0 O=a+b;                 PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   xor2    1.0 O=a^b;                 PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   or3     1.0 O=a+b+c;               PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   maj3    1.0 O=(a*b)+(b*c)+(a*c);   PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0";

struct window_manager_params : rinox::windowing::default_window_manager_params {
  static constexpr uint32_t max_num_leaves = 16u;
};

TEST_CASE("Rewiring analysis for reconvergent network",
          "[rewire_dependencies]") {
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
  fs.push_back(signal{0, 0});     // 0
  fs.push_back(signal{1, 0});     // 1
  fs.push_back(ntk.create_pi());  // 2
  fs.push_back(ntk.create_pi());  // 3
  fs.push_back(ntk.create_pi());  // 4
  fs.push_back(ntk.create_pi());  // 5
  fs.push_back(ntk.create_pi());  // 6
  fs.push_back(ntk.create_pi());  // 7

  fs.push_back(ntk.create_node(std::vector<signal>{fs[2], fs[3]}, 0));  // 8
  fs.push_back(ntk.create_node(std::vector<signal>{fs[2], fs[3]}, 1));  // 9
  fs.push_back(ntk.create_node(std::vector<signal>{fs[2], fs[3]}, 2));  // 10
  fs.push_back(ntk.create_node(std::vector<signal>{fs[4], fs[5]}, 0));  // 11
  fs.push_back(ntk.create_node(std::vector<signal>{fs[4], fs[5]}, 1));  // 12
  fs.push_back(ntk.create_node(std::vector<signal>{fs[4], fs[5]}, 2));  // 13
  fs.push_back(ntk.create_node(std::vector<signal>{fs[6], fs[7]}, 0));  // 14
  fs.push_back(ntk.create_node(std::vector<signal>{fs[6], fs[7]}, 1));  // 15
  fs.push_back(ntk.create_node(std::vector<signal>{fs[6], fs[7]}, 2));  // 16
  fs.push_back(
      ntk.create_node(std::vector<signal>{fs[9], fs[12], fs[15]}, 4));  // 17
  fs.push_back(
      ntk.create_node(std::vector<signal>{fs[8], fs[17], fs[11]}, 3));    // 18
  fs.push_back(ntk.create_node(std::vector<signal>{fs[18], fs[14]}, 1));  // 19

  ntk.create_po(fs[19]);
  ntk.create_po(fs[10]);
  ntk.create_po(fs[13]);
  ntk.create_po(fs[16]);

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk(ntk);

  window_manager_params ps;
  ps.max_num_levels_tfo = 4u;

  rinox::windowing::window_manager<DNtk, window_manager_params> window(dntk, ps,
                                                                       st);
  CHECK(window.run(dntk.get_node(fs[17])));
  rinox::windowing::window_simulator sim(dntk);
  sim.run(window);
  auto tta = sim.get(fs[2]);
  auto ttb = sim.get(fs[3]);
  auto ttc = sim.get(fs[4]);
  auto ttd = sim.get(fs[5]);
  auto tte = sim.get(fs[6]);
  auto ttf = sim.get(fs[7]);
  auto care = sim.get_careset();
  auto const exp = ~((tta & ttb) | (ttc & ttd) | (tte & ttf));
  CHECK(kitty::equal(exp, care));

  rinox::dependency::rewire_dependencies dep(ntk);
  dep.run(window, sim);

  dep.foreach_cut([&](auto const& cut, auto i) {
    for (auto l : cut.leaves) std::cout << l << " ";
    std::cout << std::endl;
  });
}
