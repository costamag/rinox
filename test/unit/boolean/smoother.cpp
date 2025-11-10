#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <kitty/kitty.hpp>
#include <rinox/boolean/metrics.hpp>
#include <rinox/boolean/smoother.hpp>
#include <rinox/io/csv/csv.hpp>
#include <rinox/math/math.hpp>
#include <rinox/synthesis/lut_decomposer.hpp>

using namespace rinox::synthesis;

#include <string>
#include <vector>

TEST_CASE("Synthesis of random partially specified Boolean functions",
          "[smoother]") {
  using CSTT = kitty::static_truth_table<7u>;
  using ISTT = kitty::ternary_truth_table<CSTT>;

  rinox::boolean::Smoother<CSTT, uint32_t, 7u> synthesizer;

  std::array<CSTT, 7u> xs;
  for (uint32_t i = 0; i < 7u; ++i) kitty::create_nth_var(xs[i], i);

  rinox::math::static_vector<CSTT const*, 7u> tts;
  rinox::math::static_vector<uint32_t, 7u> vars;
  for (uint32_t i = 0; i < 7u; ++i) {
    tts.push_back(&xs[i]);
    vars.push_back(i);
  }

  double smoothness_fast = 0;
  double smoothness_cart = 0;
  double smoothness_c4d5 = 0;
  double smoothness_prop = 0;

  std::vector<double> data_fast;
  std::vector<double> data_cart;
  std::vector<double> data_c4d5;
  std::vector<double> data_prop;

  uint32_t num_tests = 100;

  for (int i = 0; i < num_tests; ++i) {
    CSTT bits, care;
    kitty::create_random(bits);
    kitty::create_random(care);
    ISTT func(bits, care);
    auto const bdd_func_fast = synthesizer.RunFast(func, vars, tts);
    auto const bdd_func_cart = synthesizer.RunCART(func, vars, tts);
    auto const bdd_func_c4d5 = synthesizer.RunC4d5(func, vars, tts);
    auto const bdd_func_prop = synthesizer.RunLabelPropagation(func);
    CHECK(kitty::equal(bits & care, bdd_func_fast & care));
    CHECK(kitty::equal(bits & care, bdd_func_cart & care));
    CHECK(kitty::equal(bits & care, bdd_func_c4d5 & care));
    CHECK(kitty::equal(bits & care, bdd_func_prop & care));

    data_fast.push_back(rinox::boolean::Smoothness(bdd_func_fast));
    data_cart.push_back(rinox::boolean::Smoothness(bdd_func_cart));
    data_c4d5.push_back(rinox::boolean::Smoothness(bdd_func_c4d5));
    data_prop.push_back(rinox::boolean::Smoothness(bdd_func_prop));
    smoothness_fast +=
        rinox::boolean::Smoothness(bdd_func_fast) / ((double)num_tests);
    smoothness_cart +=
        rinox::boolean::Smoothness(bdd_func_cart) / ((double)num_tests);
    smoothness_c4d5 +=
        rinox::boolean::Smoothness(bdd_func_c4d5) / ((double)num_tests);
    smoothness_prop +=
        rinox::boolean::Smoothness(bdd_func_prop) / ((double)num_tests);
  }
  CHECK(smoothness_fast < smoothness_cart);
  CHECK(smoothness_cart < smoothness_c4d5);
  CHECK(smoothness_c4d5 < smoothness_prop);

  rinox::io::csv::WriteCsv("bdd_smoothness.csv", {{"fast", data_fast},
                                                  {"cart", data_cart},
                                                  {"c4d5", data_c4d5},
                                                  {"prop", data_prop}});
}
