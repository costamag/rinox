#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <rinox/network/network.hpp>
#include <rinox/opto/cell/rewrite/cell_rewrite_area.hpp>

std::string const test_library =
    "GATE   and2    1.0 O=a*b;                 PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   or2     1.0 O=a+b;                 PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   xor2    0.5 O=a^b;                 PIN * INV 1   999 3.0 0.0 3.0 "
    "0.0\n"
    "GATE   or3     1.0 O=a+b+c;               PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   and3    1.0 O=((a*b)*c);           PIN * INV 1   999 3.0 0.0 3.0 "
    "0.0\n"
    "GATE   maj3    1.0 O=(a*b)+(b*c)+(a*c);   PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   fa      1.0 C=a*b+a*c+b*c;         PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   fa      1.0 S=a^b^c;               PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   nand2   1.0 O=!(a*b);              PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   inv1    1.0 O=!a;                  PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   cxg1    1.0 O=(((!a)*b)*c);        PIN * INV 1   999 3.0 0.0 3.0 "
    "0.0";

struct CellRewriteAreaParams1 : rinox::opto::cell::RewriteAreaParamsDefault {
  static constexpr uint32_t kMaxNumLeaves = 8;
  static constexpr uint32_t num_vars_sign = 8;
  static constexpr uint32_t kInpsCapacity = 512u;
  static constexpr uint32_t kDivsCapacity = 2048u;
  static constexpr uint32_t kMffcCapacity = 512u;
  static constexpr uint32_t kTfosCapacity = 512u;
  static constexpr uint32_t kOutsCapacity = 512u;
  static constexpr uint32_t kMaxCutsSize = 6u;
  static constexpr uint32_t max_cuts_size = 6u;
  static constexpr uint32_t kMaxCutReconv = 3u;
  static constexpr uint32_t max_reconv = 3u;
  static constexpr uint32_t kMaxNumVarsDb = 6u;

  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 10;
  int max_num_levels_tfi = 10;
  int max_num_levels_tfo = 0;
};

TEST_CASE("Cell Rewrite Area - One single-output node, delay unconstrained",
          "[cell_rewrite_area]") {
  using Ntk =
      rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED,
                                    2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in(test_library);
  auto result = lorina::read_genlib(in, mockturtle::genlib_reader(gates));
  CHECK(result == lorina::return_code::success);

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED>
      lib(gates);

  static constexpr uint32_t MaxNumVars = 6u;
  using Db = rinox::databases::mapped_database<Ntk, MaxNumVars>;
  Db db(lib);

  for (int i = 0; i < 3; ++i) {
    Ntk ntk(gates);
    auto const a = ntk.create_pi();
    auto const b = ntk.create_pi();
    auto const c = ntk.create_pi();
    signal f2 = a;
    signal f3 = b;
    signal f4 = c;
    if (i == 1) std::swap(f2, f3);
    if (i == 2) std::swap(f3, f4);
    auto const f5 = ntk.create_node({f3}, 9u);
    auto const f6 = ntk.create_node({f2, f4}, 0u);
    auto const f7 = ntk.create_node({f5, f6}, 0u);
    ntk.create_po(f7);

    using DNtk = mockturtle::depth_view<Ntk>;
    rinox::windowing::window_manager_stats st;
    DNtk dntk(ntk);

    CHECK(ntk.area() == 3);
    CellRewriteAreaParams1 ps;
    ps.max_num_levels_tfi = 10;
    ps.max_num_levels_tfo = 0;
    rinox::opto::cell::RewriteArea<Ntk, Db, CellRewriteAreaParams1>(ntk, db,
                                                                    ps);
    CHECK(ntk.area() == 1);

    rinox::analyzers::evaluators::TimesEvaluator times(ntk);
    times.Run();
    CHECK(times.GetWorstArrival() == 3);
  }
}

std::string const test_library2 =
    "GATE   and2    1.0 O=a*b;                 PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   xor2    0.0 O=a^b;                 PIN * INV 1   999 3.0 0.0 3.0 "
    "0.0\n"
    "GATE   and3    1.0 O=((a*b)*c);           PIN * INV 1   999 3.0 0.0 3.0 "
    "0.0\n"
    "GATE   inv1    1.0 O=!a;                  PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   cxg1    2.0 O=(((!a)*b)*c);        PIN * INV 1   999 3.0 0.0 3.0 "
    "0.0";

TEST_CASE("Cell Rewrite Area - One single-output node, delay constrained",
          "[cell_rewrite_area]") {
  using Ntk =
      rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED,
                                    2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in(test_library2);
  auto result = lorina::read_genlib(in, mockturtle::genlib_reader(gates));
  CHECK(result == lorina::return_code::success);

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED>
      lib(gates);

  static constexpr uint32_t MaxNumVars = 6u;
  using Db = rinox::databases::mapped_database<Ntk, MaxNumVars>;
  Db db(lib);

  for (int i = 0; i < 3; ++i) {
    Ntk ntk(gates);
    auto const a = ntk.create_pi();
    auto const b = ntk.create_pi();
    auto const c = ntk.create_pi();
    signal f2 = a;
    signal f3 = b;
    signal f4 = c;
    if (i == 1) std::swap(f2, f3);
    if (i == 2) std::swap(f3, f4);
    auto const f5 = ntk.create_pi();
    auto const f6 = ntk.create_pi();
    auto const f7 = ntk.create_node({f3}, 3u);
    auto const f8 = ntk.create_node({f2, f4}, 0u);
    auto const f9 = ntk.create_node({f7, f8}, 0u);
    auto const f10 = ntk.create_node({f5, f6}, 1u);
    auto const f11 = ntk.create_node({f10, f9}, 1u);
    ntk.create_po(f11);

    using DNtk = mockturtle::depth_view<Ntk>;
    rinox::windowing::window_manager_stats st;
    DNtk dntk(ntk);

    CHECK(ntk.area() == 3);
    CellRewriteAreaParams1 ps;
    ps.max_num_levels_tfi = 10;
    ps.max_num_levels_tfo = 0;
    ps.delay_constrained = true;
    rinox::opto::cell::RewriteArea<Ntk, Db, CellRewriteAreaParams1>(ntk, db,
                                                                    ps);
    CHECK(ntk.area() == 2);

    rinox::analyzers::evaluators::TimesEvaluator times(ntk);
    times.Run();
    CHECK(times.GetWorstArrival() == 6);
  }
}