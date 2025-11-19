#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <rinox/network/network.hpp>
#include <rinox/opto/array/window_resub/array_window_resub_area.hpp>

struct Params
    : rinox::opto::array::WindowResubAreaParamsDefault {
  static constexpr uint32_t kMaxNumLeaves = RINOX_MAX_NUM_LEAVES;
  static constexpr uint32_t kMaxCutsSize = RINOX_MAX_CUTS_SIZE;
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr uint32_t kMaxInpsSolo = 6u;
  uint32_t stop_after_nth = std::numeric_limits<uint32_t>::max();
  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 10;
  int32_t odc_levels = 0;
  bool try_0resub = false;
  bool try_1resub = false;
  bool try_2resub = false;
};

TEST_CASE( "Array Window Resub | Resynthesis via squeezing", "[array_window_resub]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 2>;
  using signal = typename Ntk::signal;

  static constexpr uint32_t MaxNumVars = 6u;
  Ntk ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();
  kitty::dynamic_truth_table xor2( 2u );
  kitty::create_from_binary_string( xor2, "0110" );
  auto const f1 = ntk.create_node( { a, b }, xor2 );
  auto const f2 = ntk.create_node( { f1, c }, xor2 );
  auto const f3 = ntk.create_node( { a, f2 }, xor2 );
  auto const f4 = ntk.create_node( { b, c }, xor2 );
  ntk.create_po( f3 );
  ntk.create_po( f4 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );
  Params ps;

  ps.stop_after_nth = std::numeric_limits<uint32_t>::max();
  ps.max_num_divisors = 1024;
  ps.max_reconvergence = 10;
  ps.odc_levels = 3;
  ps.try_0resub = true;
  ps.try_1resub = false;
  ps.try_2resub = false;

  rinox::opto::array::WindowResubAreaStatsDefault st;
  rinox::opto::array::WindowResubArea<Ntk, Params>(ntk, ps, &st);
  CHECK( ntk.area() == 1 );
}

TEST_CASE( "Array Window Resub | 0-resub of a single signal", "[array_window_resub]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 2>;
  using signal = typename Ntk::signal;

  static constexpr uint32_t MaxNumVars = 6u;
  Ntk ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();
  kitty::dynamic_truth_table xor2( 2u );
  kitty::dynamic_truth_table and2( 2u );
  kitty::create_from_binary_string( xor2, "0110" );
  kitty::create_from_binary_string( and2, "1000" );
  auto const f1 = ntk.create_node( { a, b }, xor2 );
  auto const f2 = ntk.create_node( { f1, c }, xor2 );
  auto const f3 = ntk.create_node( { a, f2 }, {xor2, and2} );
  auto const f4 = ntk.create_node( { b, c }, xor2 );
  auto const f5 = ntk.create_node( { a, f2 }, and2 );
  ntk.create_po( f3 );
  ntk.create_po( f4 );
  ntk.create_po( f5 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );
  Params ps;

  ps.stop_after_nth = std::numeric_limits<uint32_t>::max();
  ps.max_num_divisors = 1024;
  ps.max_reconvergence = 10;
  ps.odc_levels = 0;
  ps.try_0resub = true;
  ps.try_1resub = false;
  ps.try_2resub = false;

  rinox::opto::array::WindowResubAreaStatsDefault st;
  rinox::opto::array::WindowResubArea<Ntk, Params>(ntk, ps, &st);
  CHECK( ntk.area() == 3 );
}
