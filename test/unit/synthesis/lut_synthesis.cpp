#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <kitty/kitty.hpp>

#include <rinox/synthesis/lut_decomposer.hpp>
#include <rinox/synthesis/synthesis.hpp>

using namespace rinox::synthesis;
using namespace rinox::evaluation;
using namespace rinox::evaluation::chains;
using namespace mockturtle;

template<uint32_t NumVars, uint32_t CutSize, typename Func, bool ExacSuppMin, typename Times>
void test_lut_dec( Func const& func, Times const& times, std::vector<std::vector<uint8_t>> const& supps, std::vector<Func> const& funcs, bool try_spfd = false, bool spfd_exact = false )
{
  lut_decomposer_params ps;
  ps.try_spfd_decompose = try_spfd;
  ps.exact_spfd = spfd_exact;
  lut_decomposer<CutSize, NumVars, ExacSuppMin> decomposer( ps );
  CHECK( decomposer.run( func, times ) );

  int j = 0;
  bool res = decomposer.foreach_spec( [&]( auto const specs, auto i ) {
    {
      std::vector<uint8_t> inputs;
      for ( auto k : specs[i].inputs )
        inputs.push_back( k );
      CHECK( inputs == supps[j] );
      CHECK( kitty::equal( specs[i].sim, funcs[j++] ) );
      return true;
    }
  } );
  CHECK( res );
}

#if 0
TEST_CASE( "LUT decomposer for parity 4 3- shannon", "[lut_synthesis]" )
{
  static constexpr uint32_t num_vars = 3u;
  static constexpr uint32_t cut_size = 4u;

  using ctt_s = kitty::static_truth_table<num_vars>;
  using itt_s = kitty::ternary_truth_table<ctt_s>;
  using ctt_l = kitty::static_truth_table<cut_size>;
  using itt_l = kitty::ternary_truth_table<ctt_l>;

  ctt_l ctt;
  kitty::create_parity( ctt );
  itt_l func( ctt );
  std::vector<double> times{ 0, 0, 0, 3 };

  std::vector<std::vector<uint8_t>> supps;
  supps = { { 0, 1, 2 }, { 0, 1, 2 }, { 3, 4, 5 } };
  std::vector<ctt_l> ctts_l( 3u );
  kitty::create_from_binary_string( ctts_l[0], "1001011010010110" );
  kitty::create_from_binary_string( ctts_l[1], "0110100101101001" );
  kitty::create_from_binary_string( ctts_l[2], "0110100110010110" );
  std::vector<itt_l> funcs;
  funcs.emplace_back( ctts_l[0] );
  funcs.emplace_back( ctts_l[1] );
  funcs.emplace_back( ctts_l[2] );
  test_lut_dec<num_vars, cut_size, itt_l, /* exact support minimization */ true>( func, times, supps, funcs );
  test_lut_dec<num_vars, cut_size, itt_l, /* exact support minimization */ false>( func, times, supps, funcs );
}

TEST_CASE( "LUT decomposer for parity 5 3- shannon", "[lut_synthesis]" )
{
  static constexpr uint32_t num_vars = 3u;
  static constexpr uint32_t cut_size = 5u;

  using ctt_s = kitty::static_truth_table<num_vars>;
  using itt_s = kitty::ternary_truth_table<ctt_s>;
  using ctt_l = kitty::static_truth_table<cut_size>;
  using itt_l = kitty::ternary_truth_table<ctt_l>;

  ctt_l ctt;
  kitty::create_parity( ctt );
  itt_l func( ctt );
  std::vector<double> times{ 0, 0, 0, 0, 3 };

  std::vector<std::vector<uint8_t>> supps;
  supps = { { 0, 1, 2 }, { 0, 1, 2 }, { 3, 5, 6 }, { 0, 1, 2 }, { 0, 1, 2 }, { 3, 8, 9 }, { 4, 7, 10 } };
  std::vector<ctt_l> ctts_l( 7u );
  kitty::create_from_binary_string( ctts_l[0], "10010110100101101001011010010110" );
  kitty::create_from_binary_string( ctts_l[1], "01101001011010010110100101101001" );
  kitty::create_from_binary_string( ctts_l[2], "01101001100101100110100110010110" );
  kitty::create_from_binary_string( ctts_l[3], "01101001011010010110100101101001" );
  kitty::create_from_binary_string( ctts_l[4], "10010110100101101001011010010110" );
  kitty::create_from_binary_string( ctts_l[5], "10010110011010011001011001101001" );
  kitty::create_from_binary_string( ctts_l[6], "10010110011010010110100110010110" );
  std::vector<itt_l> funcs;
  funcs.emplace_back( ctts_l[0] );
  funcs.emplace_back( ctts_l[1] );
  funcs.emplace_back( ctts_l[2] );
  funcs.emplace_back( ctts_l[3] );
  funcs.emplace_back( ctts_l[4] );
  funcs.emplace_back( ctts_l[5] );
  funcs.emplace_back( ctts_l[6] );
  test_lut_dec<num_vars, cut_size, itt_l, /* exact support minimization */ true>( func, times, supps, funcs );
  test_lut_dec<num_vars, cut_size, itt_l, /* exact support minimization */ false>( func, times, supps, funcs );
}

TEST_CASE( "LUT decomposer for parity 5 3- SPFD enumeration", "[lut_synthesis]" )
{
  static constexpr uint32_t num_vars = 3u;
  static constexpr uint32_t cut_size = 5u;

  using ctt_s = kitty::static_truth_table<num_vars>;
  using itt_s = kitty::ternary_truth_table<ctt_s>;
  using ctt_l = kitty::static_truth_table<cut_size>;
  using itt_l = kitty::ternary_truth_table<ctt_l>;

  ctt_l ctt;
  kitty::create_parity( ctt );
  itt_l func( ctt );
  std::vector<double> times{ 0, 0, 0, 0, 3 };

  std::vector<std::vector<uint8_t>> supps;
  supps = { { 0, 1, 2 }, { 3, 4, 5 } };
  std::vector<ctt_l> ctts_l( 7u );
  kitty::create_from_binary_string( ctts_l[0], "10010110100101101001011010010110" );
  kitty::create_from_binary_string( ctts_l[1], "10010110011010010110100110010110" );
  std::vector<itt_l> funcs;
  funcs.emplace_back( ctts_l[0] );
  funcs.emplace_back( ctts_l[1] );
  test_lut_dec<num_vars, cut_size, itt_l, false>( func, times, supps, funcs, true, false );
  test_lut_dec<num_vars, cut_size, itt_l, true>( func, times, supps, funcs, true, false );
  test_lut_dec<num_vars, cut_size, itt_l, true>( func, times, supps, funcs, true, true );
  test_lut_dec<num_vars, cut_size, itt_l, false>( func, times, supps, funcs, true, true );
}
#endif