#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <kitty/kitty.hpp>
#include <rinox/math/math.hpp>
#include <rinox/network/network.hpp>
#include <rinox/synthesis/lut/spin_dec.hpp>

struct LutSynthParams1 : rinox::synthesis::lut::DefaultLutSynthParams {
  static constexpr uint32_t kMaxCutsSize = 5u;
  static constexpr uint32_t kMaxInpsSolo = 4u;
  static constexpr uint32_t kMaxInpsDual = 3u;
  static constexpr uint32_t kMaxNumOutputs = 1u;
  static constexpr uint32_t kMaxNumGates = 32u;
  static constexpr bool kExactSuppMin = false;
};

struct LutSynthParams2 : rinox::synthesis::lut::DefaultLutSynthParams {
  static constexpr uint32_t kMaxCutsSize = 5u;
  static constexpr uint32_t kMaxInpsSolo = 4u;
  static constexpr uint32_t kMaxInpsDual = 3u;
  static constexpr uint32_t kMaxNumOutputs = 2u;
  static constexpr uint32_t kMaxNumGates = 32u;
  static constexpr bool kExactSuppMin = false;
};

TEST_CASE("Spin-based decomposition - 1 single-output LUT", "[lut_spin_dec]") {
  rinox::synthesis::lut::SpinDecomposition<LutSynthParams1> synth;
  using CSFunc = kitty::static_truth_table<LutSynthParams1::kMaxCutsSize>;
  using ISFunc = kitty::ternary_truth_table<CSFunc>;
  CSFunc csfunc;
  std::array<CSFunc, LutSynthParams1::kMaxCutsSize> xs;
  for (int i = 0; i < LutSynthParams1::kMaxCutsSize; ++i) {
    kitty::create_nth_var(xs[i], i);
  }
  ISFunc func;
  // f = x0 & x1 & x2
  func._care = xs[0] | ~xs[0];
  func._bits = xs[0] & xs[1] & xs[2] & xs[4];
  bool res = synth.Run(func);
  CHECK(res);

  using Ntk =
      rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED,
                                    1>;
  Ntk ntk;
  std::vector<typename Ntk::signal> pis;
  pis.push_back(ntk.create_pi());
  pis.push_back(ntk.create_pi());
  pis.push_back(ntk.create_pi());
  pis.push_back(ntk.create_pi());
  pis.push_back(ntk.create_pi());
  auto const sigs = synth.Insert(ntk, pis);
  for (auto const s : sigs) {
    std::cout << s << std::endl;
  }
}

TEST_CASE("Spin-based decomposition - 1 dual-output LUT", "[lut_spin_dec]") {
  rinox::synthesis::lut::SpinDecomposition<LutSynthParams2> synth;
  using CSFunc = kitty::static_truth_table<LutSynthParams2::kMaxCutsSize>;
  using ISFunc = kitty::ternary_truth_table<CSFunc>;
  CSFunc csfunc;
  std::array<CSFunc, LutSynthParams2::kMaxCutsSize> xs;
  for (int i = 0; i < LutSynthParams2::kMaxCutsSize; ++i) {
    kitty::create_nth_var(xs[i], i);
  }
  ISFunc func1, func2;
  // f = x0 & x1 & x2
  func1._care = xs[0] | ~xs[0];
  func1._bits = xs[0] & xs[1] & xs[2] & xs[4];

  func2._care = xs[0] | ~xs[0];
  func2._bits = xs[0] | xs[1] | xs[2];

  rinox::math::static_vector<ISFunc, LutSynthParams2::kMaxNumOutputs> funcs;
  funcs.push_back(func1);
  funcs.push_back(func2);
  bool res = synth.Run(funcs);
  CHECK(res);
}