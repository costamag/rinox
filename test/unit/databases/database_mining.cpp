#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <lorina/genlib.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <rinox/analyzers/trackers/arrival_times_tracker.hpp>
#include <rinox/databases/database_mining.hpp>
#include <rinox/databases/mapped_database.hpp>
#include <rinox/evaluation/chains/bound_chain.hpp>
#include <rinox/network/network.hpp>
#include <sstream>
#include <vector>

using namespace mockturtle;

std::string const test_library =
    "GATE   and2    1.0 O=a*b;                 PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0\n"
    "GATE   xor2    0.0 O=a^b;                 PIN * INV 1   999 3.0 0.0 3.0 "
    "0.0\n"
    "GATE   inv1    1.0 O=!a;                  PIN * INV 1   999 1.0 0.0 1.0 "
    "0.0";

TEST_CASE("Mining from AIG - Basic addition", "[database_mining]") {
  using bound_network =
      rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED,
                                    2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in(test_library);
  auto result = lorina::read_genlib(in, genlib_reader(gates));
  CHECK(result == lorina::return_code::success);

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED>
      lib(gates);
  static constexpr uint32_t MaxNumVars = 6u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db(lib);

  rinox::databases::MiningParams ps;
  rinox::databases::MiningStats st;
  rinox::databases::DatabaseMining<bound_network> miner(gates, ps, st);

  bound_network ntk(gates);
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();

  auto const f1 = ntk.create_node({a, b}, 0u);    // and2
  auto const f2 = ntk.create_node({c, d}, 1u);    // xor2
  auto const f3 = ntk.create_node({f1, f2}, 0u);  // and2
  auto const f4 = ntk.create_node({f3, e}, 1u);   // xor2
  ntk.create_po(f4);

  miner.Collect(ntk);
  CHECK(st.num_added == 7);
  miner.Collect(ntk);
  CHECK(st.num_added == 7);

  auto const f5 = ntk.create_node({c, b}, 0u);    // and2
  auto const f6 = ntk.create_node({a, d}, 1u);    // xor2
  auto const f7 = ntk.create_node({f5, f6}, 0u);  // and2
  ntk.create_po(f7);

  miner.Collect(ntk);
  CHECK(st.num_added == 7);

  auto const f8 = ntk.create_node({a, b}, 1u);     // xor2
  auto const f9 = ntk.create_node({c, d}, 1u);     // xor2
  auto const f10 = ntk.create_node({f8, f9}, 1u);  // xor2
  ntk.create_po(f10);

  miner.Collect(ntk);
  CHECK(st.num_added == 9);

  auto const f11 = ntk.create_node({a}, 2u);       // inv1 +0
  auto const f12 = ntk.create_node({b}, 2u);       // inv1 +0
  auto const f13 = ntk.create_node({f11, b}, 0u);  // and2 +1

  ntk.create_po(f13);
  miner.Collect(ntk);
  CHECK(st.num_added == 10);

  auto const f14 = ntk.create_node({f12, a}, 0u);  // and2 +0
  ntk.create_po(f14);
  miner.Collect(ntk);
  CHECK(st.num_added == 10);

  auto const f15 = ntk.create_node({f13}, 2u);  // inv1 +2
  ntk.create_po(f15);
  miner.Collect(ntk);
  CHECK(st.num_added == 12);

  auto const f16 = ntk.create_node({f14}, 2u);  // inv1 +0
  ntk.create_po(f16);
  miner.Collect(ntk);
  CHECK(st.num_added == 12);

  auto const f17 = ntk.create_node({f15, f16}, 0u);  // inv1 +2
  ntk.create_po(f17);
  miner.Collect(ntk);
  CHECK(st.num_added == 20);

  auto const f18 = ntk.create_node({f17}, 2u);  // inv1 +3
  ntk.create_po(f18);
  miner.Collect(ntk);
  CHECK(st.num_added == 27);
}
