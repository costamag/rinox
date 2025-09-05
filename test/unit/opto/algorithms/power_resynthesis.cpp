#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <rinox/analyzers/evaluators/power_evaluator.hpp>
#include <rinox/network/network.hpp>
#include <rinox/opto/algorithms/resynthesize.hpp>

std::string const test_library = "GATE   and2    1.0 O=a*b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 0
                                 "GATE   or2     1.0 O=a+b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 1
                                 "GATE   xor2    0.5 O=a^b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 2
                                 "GATE   or3     1.0 O=a+b+c;               PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 3
                                 "GATE   and3    1.0 O=((a*b)*c);           PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 4
                                 "GATE   maj3    1.0 O=(a*b)+(b*c)+(a*c);   PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 5
                                 "GATE   fa      1.0 C=a*b+a*c+b*c;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 6
                                 "GATE   fa      1.0 S=a^b^c;               PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 7
                                 "GATE   nand2   1.0 O=!(a*b);              PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 8
                                 "GATE   inv1    1.0 O=!a;                  PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 9
                                 "GATE   xor3    1.0 O=a^b^c;               PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 10
                                 "GATE   nor3    1.0 O=!(a+b+c);            PIN * INV 1   999 1.0 0.0 1.0 0.0";  // 11

struct custom_power_rewire_params : rinox::opto::algorithms::default_resynthesis_params<8u>
{
  bool try_rewire = true;
  bool try_struct = false;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  int32_t odc_levels = 0;
  static constexpr uint32_t max_cuts_size = 6u;
};

TEST_CASE( "Power resynthesis via rewiring - single-output gate without don't cares", "[power_resynthesis]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  using Db = rinox::databases::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const f1 = ntk.create_node( { a }, 9u );
  auto const f2 = ntk.create_node( { f1, b }, 2u );
  auto const f3 = ntk.create_node( { f2 }, 9u );
  auto const f4 = ntk.create_node( { f3, c, d }, 5u );
  auto const f5 = ntk.create_node( { a, b }, 2u );
  auto const f6 = ntk.create_node( { f5, c }, 0u );

  ntk.create_po( f4 );
  ntk.create_po( f6 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );

  std::vector<kitty::static_truth_table<8>> tts_ini( 4u );
  std::vector<kitty::static_truth_table<8>> tts_end( 4u );
  auto i = 0u;
  for ( ; i < 4; ++i )
    kitty::create_nth_var( tts_ini[i], i );
  for ( ; i < 8; ++i )
    kitty::create_nth_var( tts_end[i - 4], i );
  rinox::analyzers::utils::workload<kitty::static_truth_table<8>, 10> workload( tts_ini, tts_end );
  rinox::analyzers::evaluators::power_evaluator_stats st;
  rinox::analyzers::evaluators::power_evaluator<DNtk, kitty::static_truth_table<8>, 10> power( dntk, st );

  power.run( workload );
  auto const power_before = st.dyn_power;
  custom_power_rewire_params ps;
  rinox::opto::algorithms::power_resynthesize<DNtk, Db, custom_power_rewire_params>( dntk, db, ps );
  power.run( workload );
  auto const power_after = st.dyn_power;
  CHECK( dntk.area() == 2.5 );

  CHECK( power_after < power_before );
}

TEST_CASE( "Power resynthesis via rewiring - single-output gate with don't cares", "[power_resynthesis]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  using Db = rinox::databases::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const f1 = ntk.create_node( { b }, 9u );
  auto const f2 = ntk.create_node( { a, f1 }, 2u );
  auto const f3 = ntk.create_node( { f2 }, 9u );
  auto const f4 = ntk.create_node( { f3, c, d }, 10u );
  auto const f5 = ntk.create_node( { a, b }, 8u );
  auto const f6 = ntk.create_node( { f4, f5 }, 0u );
  auto const f7 = ntk.create_node( { a, b }, 1u );
  auto const f8 = ntk.create_node( { f7, c }, 2u );

  ntk.create_po( f6 );
  ntk.create_po( f8 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );
  std::vector<kitty::static_truth_table<8>> tts_ini( 4u );
  std::vector<kitty::static_truth_table<8>> tts_end( 4u );
  auto i = 0u;
  for ( ; i < 4; ++i )
    kitty::create_nth_var( tts_ini[i], i );
  for ( ; i < 8; ++i )
    kitty::create_nth_var( tts_end[i - 4], i );
  rinox::analyzers::utils::workload<kitty::static_truth_table<8>, 10> workload( tts_ini, tts_end );
  rinox::analyzers::evaluators::power_evaluator_stats st;
  rinox::analyzers::evaluators::power_evaluator<DNtk, kitty::static_truth_table<8>, 10> power( dntk, st );
  power.run( workload );

  auto const power_before = st.dyn_power;
  custom_power_rewire_params ps;
  ps.window_manager_ps.odc_levels = 3;
  rinox::opto::algorithms::power_resynthesize<DNtk, Db, custom_power_rewire_params>( dntk, db, ps );

  power.run( workload );
  auto const power_after = st.dyn_power;
  CHECK( dntk.area() == 4.5 );

  CHECK( power_after < power_before );
}

TEST_CASE( "Power resynthesis via rewiring - multiple-output gate without don't cares", "[power_resynthesis]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  using Db = rinox::databases::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();
  auto const f1 = ntk.create_node( { c, d, e }, 10u );
  auto const f2 = ntk.create_node( { a, b }, 0u );
  auto const f3 = ntk.create_node( { c, d }, 2u );
  auto const f4 = ntk.create_node( { e, f3 }, 2u );
  auto const f5 = ntk.create_node( { c, f2, f4 }, { 6u, 7u } );

  ntk.create_po( f1 );
  ntk.create_po( { f5.index, 0 } );
  ntk.create_po( { f5.index, 1 } );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );
  std::vector<kitty::static_truth_table<10>> tts_ini( 5u );
  std::vector<kitty::static_truth_table<10>> tts_end( 5u );
  auto i = 0u;
  for ( ; i < 5; ++i )
    kitty::create_nth_var( tts_ini[i], i );
  for ( ; i < 10; ++i )
    kitty::create_nth_var( tts_end[i - 5], i );
  rinox::analyzers::utils::workload<kitty::static_truth_table<10>, 10> workload( tts_ini, tts_end );
  rinox::analyzers::evaluators::power_evaluator_stats st;
  rinox::analyzers::evaluators::power_evaluator<DNtk, kitty::static_truth_table<10>, 10> power( dntk, st );

  power.run( workload );
  auto const power_before = st.dyn_power;
  custom_power_rewire_params ps;
  rinox::opto::algorithms::power_resynthesize<DNtk, Db, custom_power_rewire_params>( dntk, db, ps );

  power.run( workload );
  auto const power_after = st.dyn_power;
  CHECK( dntk.area() == 3.0 );

  CHECK( power_after < power_before );
}

TEST_CASE( "Power resynthesis via rewiring - multiple-output gate with don't cares", "[power_resynthesis]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  using Db = rinox::databases::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const f1 = ntk.create_node( { a }, 9u );
  auto const f2 = ntk.create_node( { f1, b }, 2u );
  auto const f3 = ntk.create_node( { f2 }, 9u );
  auto const f4 = ntk.create_node( { f3, c, d }, { 6u, 7u } );
  auto const f5 = ntk.create_node( { a, b }, 0u );
  auto const f6 = ntk.create_node( { { f4.index, 0 }, { f4.index, 1 } }, 2u );
  auto const f7 = ntk.create_node( { f6, f5 }, 1u );
  auto const f8 = ntk.create_node( { a, b }, 1u );
  auto const f9 = ntk.create_node( { c, d }, 2u );
  auto const f10 = ntk.create_node( { f8, f9 }, 2u );

  ntk.create_po( f7 );
  ntk.create_po( f10 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );
  std::vector<kitty::static_truth_table<8>> tts_ini( 4u );
  std::vector<kitty::static_truth_table<8>> tts_end( 4u );
  auto i = 0u;
  for ( ; i < 4; ++i )
    kitty::create_nth_var( tts_ini[i], i );
  for ( ; i < 8; ++i )
    kitty::create_nth_var( tts_end[i - 4], i );
  rinox::analyzers::utils::workload<kitty::static_truth_table<8>, 10> workload( tts_ini, tts_end );
  rinox::analyzers::evaluators::power_evaluator_stats st;
  rinox::analyzers::evaluators::power_evaluator<DNtk, kitty::static_truth_table<8>, 10> power( dntk, st );

  power.run( workload );
  auto const power_before = st.dyn_power;
  custom_power_rewire_params ps;
  ps.window_manager_ps.odc_levels = 3;
  rinox::opto::algorithms::power_resynthesize<DNtk, Db, custom_power_rewire_params>( dntk, db, ps );

  power.run( workload );
  auto const power_after = st.dyn_power;
  CHECK( dntk.area() == 4.0 );

  CHECK( power_after < power_before );
}

struct custom_power_struct_params1 : rinox::opto::algorithms::default_resynthesis_params<8u>
{
  bool try_rewire = false;
  bool try_struct = true;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  int32_t odc_levels = 0;
  static constexpr uint32_t max_cuts_size = 6u;
};

struct custom_power_struct_params2 : rinox::opto::algorithms::default_resynthesis_params<8u>
{
  bool try_rewire = false;
  bool try_struct = true;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  int32_t odc_levels = 2;
  static constexpr uint32_t max_cuts_size = 6u;
};

TEST_CASE( "Power resynthesis via cut rewriting - single-output gate without don't cares", "[power_resynthesis]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  using Db = rinox::databases::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> list;
  list.add_inputs( MaxNumVars );
  list.add_output( list.add_gate( { list.add_gate( { 0, 1 }, 0 ) }, 9 ) );
  CHECK( db.add( list ) );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f1 = ntk.create_node( { a }, 9u );
  auto const f2 = ntk.create_node( { c, f1 }, 0u );
  auto const f3 = ntk.create_node( { c }, 9u );
  auto const f4 = ntk.create_node( { f2, f3 }, 2u );

  ntk.create_po( f4 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );
  std::vector<kitty::static_truth_table<6>> tts_ini( 3u );
  std::vector<kitty::static_truth_table<6>> tts_end( 3u );
  auto i = 0u;
  for ( ; i < 3; ++i )
    kitty::create_nth_var( tts_ini[i], i );
  for ( ; i < 6; ++i )
    kitty::create_nth_var( tts_end[i - 3], i );
  rinox::analyzers::utils::workload<kitty::static_truth_table<6>, 10> workload( tts_ini, tts_end );
  rinox::analyzers::evaluators::power_evaluator_stats st;
  rinox::analyzers::evaluators::power_evaluator<DNtk, kitty::static_truth_table<6>, 10> power( dntk, st );

  power.run( workload );
  auto const power_before = st.dyn_power;
  custom_power_struct_params1 ps;
  ps.window_manager_ps.odc_levels = 0;
  rinox::opto::algorithms::power_resynthesize<DNtk, Db, custom_power_struct_params1>( dntk, db, ps );

  power.run( workload );

  auto const power_after = st.dyn_power;
  CHECK( dntk.area() == 2 );

  CHECK( power_after < power_before );
}

struct custom_power_struct_params3 : rinox::opto::algorithms::default_resynthesis_params<6u>
{
  bool try_rewire = false;
  bool try_struct = true;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  static constexpr uint32_t max_cuts_size = 3u;
};

TEST_CASE( "Power resynthesis via cut rewriting - single-output gate with don't cares", "[power_resynthesis]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 3u;
  using Db = rinox::databases::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> list;
  list.add_inputs( MaxNumVars );
  list.add_output( list.add_gate( { 0, 1, 2 }, 5 ) );
  CHECK( db.add( list ) );

  Ntk ntk( gates );
  auto const c = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const a = ntk.create_pi();
  auto const f0 = ntk.create_node( { b }, 9u );
  auto const f1 = ntk.create_node( { a, f0 }, 2u );
  auto const f2 = ntk.create_node( { f1, c }, 2u );
  auto const f4 = ntk.create_node( { a, b, c }, 4u );
  auto const f5 = ntk.create_node( { a, b, c }, 11u );
  auto const f6 = ntk.create_node( { f4, f5 }, 1u );
  auto const f7 = ntk.create_node( { f6, f2 }, 1u );

  ntk.create_po( f7 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );
  std::vector<kitty::static_truth_table<6>> tts_ini( 3u );
  std::vector<kitty::static_truth_table<6>> tts_end( 3u );
  auto i = 0u;
  for ( ; i < 3; ++i )
    kitty::create_nth_var( tts_ini[i], i );
  for ( ; i < 6; ++i )
    kitty::create_nth_var( tts_end[i - 3], i );
  rinox::analyzers::utils::workload<kitty::static_truth_table<6>, 10> workload( tts_ini, tts_end );
  rinox::analyzers::evaluators::power_evaluator_stats st;
  rinox::analyzers::evaluators::power_evaluator<DNtk, kitty::static_truth_table<6>, 10> power( dntk, st );

  power.run( workload );
  auto const power_before = st.dyn_power;
  custom_power_struct_params3 ps;
  ps.window_manager_ps.odc_levels = 3;
  rinox::opto::algorithms::power_resynthesize<DNtk, Db, custom_power_struct_params3>( dntk, db, ps );

  power.run( workload );

  auto const power_after = st.dyn_power;
  CHECK( dntk.area() == 5 );

  CHECK( power_after < power_before );
}

struct custom_power_window_params1 : rinox::opto::algorithms::default_resynthesis_params<6u>
{
  bool try_rewire = false;
  bool try_struct = false;
  bool try_window = true;
  bool try_simula = false;
  bool dynamic_database = false;
  static constexpr uint32_t max_cuts_size = 3u;
};

TEST_CASE( "Power resynthesis via window rewriting - single-output gate without don't cares", "[power_resynthesis]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 3u;
  using Db = rinox::databases::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> list;
  list.add_inputs( MaxNumVars );
  list.add_output( list.add_gate( { list.add_gate( { 0, 1 }, 1 ), 2 }, 1 ) );
  CHECK( db.add( list ) );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f0 = ntk.create_node( { a }, 9 );
  auto const f1 = ntk.create_node( { f0, b }, 0u );
  auto const f2 = ntk.create_node( { f0, c }, 0u );
  auto const f3 = ntk.create_node( { b, c }, 0u );
  auto const f4 = ntk.create_node( { f1, f2 }, 2u );
  auto const f5 = ntk.create_node( { f4, f3 }, 2u );

  ntk.create_po( f1 );
  ntk.create_po( f2 );
  ntk.create_po( f3 );
  ntk.create_po( f5 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );
  std::vector<kitty::static_truth_table<6>> tts_ini( 3u );
  std::vector<kitty::static_truth_table<6>> tts_end( 3u );
  auto i = 0u;
  for ( ; i < 3; ++i )
    kitty::create_nth_var( tts_ini[i], i );
  for ( ; i < 6; ++i )
    kitty::create_nth_var( tts_end[i - 3], i );
  rinox::analyzers::utils::workload<kitty::static_truth_table<6>, 10> workload( tts_ini, tts_end );
  rinox::analyzers::evaluators::power_evaluator_stats st;
  rinox::analyzers::evaluators::power_evaluator<DNtk, kitty::static_truth_table<6>, 10> power( dntk, st );

  power.run( workload );
  auto const power_before = st.dyn_power;
  custom_power_window_params1 ps;
  ps.window_manager_ps.odc_levels = 0;
  rinox::opto::algorithms::power_resynthesize<DNtk, Db, custom_power_window_params1>( dntk, db, ps );

  power.run( workload );

  auto const power_after = st.dyn_power;
  CHECK( dntk.area() == 6 );

  CHECK( power_after < power_before );
}

TEST_CASE( "Power resynthesis via window rewriting - single-output gate with don't cares", "[power_resynthesis]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 3u;
  using Db = rinox::databases::mapped_database<Ntk, MaxNumVars>;
  Db db( lib );

  rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> list;
  list.add_inputs( MaxNumVars );
  list.add_output( list.add_gate( { list.add_gate( { 0 }, 9 ), 2 }, 0 ) );
  CHECK( db.add( list ) );

  Ntk ntk( gates );
  auto const a = ntk.create_pi(); // 2
  auto const b = ntk.create_pi(); // 3
  auto const c = ntk.create_pi(); // 4
  auto const d = ntk.create_pi(); // 5
  auto const f1 = ntk.create_node( { a, b }, 2 );
  auto const f2 = ntk.create_node( { c, d }, 2 );
  auto const f3 = ntk.create_node( { f1 }, 9 );
  auto const f4 = ntk.create_node( { f2, f3 }, 2 );
  auto const f5 = ntk.create_node( { f1 }, 9 );
  auto const f6 = ntk.create_node( { f2, f5 }, 0 );
  auto const f7 = ntk.create_node( { f4 }, 9 );
  auto const f8 = ntk.create_node( { f6, f7 }, 1 );

  ntk.create_po( f1 );
  ntk.create_po( f2 );
  ntk.create_po( f8 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );

  std::vector<kitty::static_truth_table<8>> tts_ini( 4u );
  std::vector<kitty::static_truth_table<8>> tts_end( 4u );
  auto i = 0u;
  for ( ; i < 4; ++i )
    kitty::create_nth_var( tts_ini[i], i );
  for ( ; i < 8; ++i )
    kitty::create_nth_var( tts_end[i - 4], i );
  rinox::analyzers::utils::workload<kitty::static_truth_table<8>, 10> workload( tts_ini, tts_end );
  rinox::analyzers::evaluators::power_evaluator_stats st;
  rinox::analyzers::evaluators::power_evaluator<DNtk, kitty::static_truth_table<8>, 10> power( dntk, st );

  power.run( workload );
  auto const power_before = st.dyn_power;
  custom_power_window_params1 ps;
  ps.window_manager_ps.odc_levels = 3;
  rinox::opto::algorithms::power_resynthesize<DNtk, Db, custom_power_window_params1>( dntk, db, ps );

  power.run( workload );

  auto const power_after = st.dyn_power;
  CHECK( dntk.area() == 6 );

  CHECK( power_after < power_before );
}