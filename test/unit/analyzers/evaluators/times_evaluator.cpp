// SPDX-License-Identifier: MIT
// Tests for rinox boolean/simd_operations.hpp

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <lorina/genlib.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <rinox/analyzers/evaluators/times_evaluator.hpp>
#include <rinox/network/network.hpp>
#include <vector>

using namespace mockturtle;
using namespace rinox;
using namespace network;
using namespace analyzers::evaluators;

std::string const test_library =
    "GATE   inv1    1 O=!a;            PIN * INV 1 999 1.0 1.0 1.0 1.0\n"
    "GATE   and2    3 O=a*b;           PIN * INV 1 999 2.0 2.0 2.0 2.0\n"
    "GATE   xor2    4 O=a^b;           PIN * UNKNOWN 2 999 3.0 3.0 3.0 3.0";

TEST_CASE("Evaluate times using tiers-based evaluator", "[times_evaluator]") {
  using bound_network = bound_network<design_type_t::CELL_BASED, 2>;
  using node = typename bound_network::node;
  using signal = typename bound_network::signal;
  std::vector<gate> gates;

  std::istringstream in(test_library);
  auto result = lorina::read_genlib(in, genlib_reader(gates));
  CHECK(result == lorina::return_code::success);

  bound_network ntk(gates);
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  /* chain of inverters */
  auto const f1 = ntk.create_node<std::vector<signal>, false>({a}, 0);
  auto const f2 = ntk.create_node<std::vector<signal>, false>({b}, 0);
  auto const f3 = ntk.create_node<std::vector<signal>, false>({c}, 0);
  auto const f4 = ntk.create_node<std::vector<signal>, false>({d}, 0);
  auto const f5 = ntk.create_node<std::vector<signal>, false>({f1, f2}, 2);
  auto const f6 = ntk.create_node<std::vector<signal>, false>({f2, f3}, 1);
  auto const f7 = ntk.create_node<std::vector<signal>, false>({f3, f4}, 1);
  auto const f8 = ntk.create_node<std::vector<signal>, false>({f6}, 0);
  auto const f9 = ntk.create_node<std::vector<signal>, false>({f5, f8}, 1);
  auto const f10 = ntk.create_node<std::vector<signal>, false>({f7, f8}, 1);

  ntk.create_po(f9);
  ntk.create_po(f10);
  TimesEvaluator<bound_network> eval(ntk);
  CHECK(eval.GetTier(ntk.get_node(f1)) == 1);
  CHECK(eval.GetTier(ntk.get_node(f2)) == 1);
  CHECK(eval.GetTier(ntk.get_node(f3)) == 1);
  CHECK(eval.GetTier(ntk.get_node(f4)) == 1);
  CHECK(eval.GetTier(ntk.get_node(f5)) == 2);
  CHECK(eval.GetTier(ntk.get_node(f6)) == 2);
  CHECK(eval.GetTier(ntk.get_node(f7)) == 2);
  CHECK(eval.GetTier(ntk.get_node(f8)) == 3);
  CHECK(eval.GetTier(ntk.get_node(f9)) == 4);
  CHECK(eval.GetTier(ntk.get_node(f10)) == 4);

  CHECK(eval.GetArrival(f1) == 1);
  CHECK(eval.GetArrival(f2) == 1);
  CHECK(eval.GetArrival(f3) == 1);
  CHECK(eval.GetArrival(f4) == 1);
  CHECK(eval.GetArrival(f5) == 4);
  CHECK(eval.GetArrival(f6) == 3);
  CHECK(eval.GetArrival(f7) == 3);
  CHECK(eval.GetArrival(f8) == 4);
  CHECK(eval.GetArrival(f9) == 6);
  CHECK(eval.GetArrival(f10) == 6);

  CHECK(eval.GetRequired(f1) == 1);
  CHECK(eval.GetRequired(f2) == 1);
  CHECK(eval.GetRequired(f3) == 1);
  CHECK(eval.GetRequired(f4) == 2);
  CHECK(eval.GetRequired(f5) == 4);
  CHECK(eval.GetRequired(f6) == 3);
  CHECK(eval.GetRequired(f7) == 4);
  CHECK(eval.GetRequired(f8) == 4);
  CHECK(eval.GetRequired(f9) == 6);
  CHECK(eval.GetRequired(f10) == 6);

  ntk.substitute_node(ntk.get_node(f8), f6);
  eval.UpdateTiers(ntk.get_node(f6), 1, 4);

  CHECK(eval.GetTier(ntk.get_node(f1)) == 1);
  CHECK(eval.GetTier(ntk.get_node(f2)) == 1);
  CHECK(eval.GetTier(ntk.get_node(f3)) == 1);
  CHECK(eval.GetTier(ntk.get_node(f4)) == 1);
  CHECK(eval.GetTier(ntk.get_node(f5)) == 2);
  CHECK(eval.GetTier(ntk.get_node(f6)) == 2);
  CHECK(eval.GetTier(ntk.get_node(f7)) == 2);
  CHECK(eval.GetTier(ntk.get_node(f8)) == 3);
  CHECK(eval.GetTier(ntk.get_node(f9)) == 3);
  CHECK(eval.GetTier(ntk.get_node(f10)) == 3);

  CHECK(eval.GetArrival(f1) == 1);
  CHECK(eval.GetArrival(f2) == 1);
  CHECK(eval.GetArrival(f3) == 1);
  CHECK(eval.GetArrival(f4) == 1);
  CHECK(eval.GetArrival(f5) == 4);
  CHECK(eval.GetArrival(f6) == 3);
  CHECK(eval.GetArrival(f7) == 3);
  CHECK(eval.GetArrival(f8) == 4);
  CHECK(eval.GetArrival(f9) == 6);
  CHECK(eval.GetArrival(f10) == 5);

  CHECK(eval.GetRequired(f1) == 1);
  CHECK(eval.GetRequired(f2) == 1);
  CHECK(eval.GetRequired(f3) == 2);
  CHECK(eval.GetRequired(f4) == 2);
  CHECK(eval.GetRequired(f5) == 4);
  CHECK(eval.GetRequired(f6) == 4);
  CHECK(eval.GetRequired(f7) == 4);
  CHECK(eval.GetRequired(f8) == 4);
  CHECK(eval.GetRequired(f9) == 6);
  CHECK(eval.GetRequired(f10) == 6);
}
