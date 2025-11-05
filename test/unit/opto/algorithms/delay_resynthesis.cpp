#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
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
                                 "GATE   inv1    1.0 O=!a;                  PIN * INV 1   999 1.0 0.0 1.0 0.0";  // 9

struct custom_delay_rewire_params : rinox::opto::algorithms::default_resynthesis_params<8u>
{
  bool try_rewire = true;
  bool try_struct = false;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  static constexpr uint32_t max_cuts_size = 6u;
  static constexpr uint32_t max_num_vars = 6u;
};

TEST_CASE( "Delay resynthesis via rewiring - single-output gate without don't cares", "[delay_resynthesis]" )
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
  auto const f1 = ntk.create_node( { a, b }, 0u );
  auto const f2 = ntk.create_node( { a, b }, 8u );
  auto const f3 = ntk.create_node( { c }, 9u );
  auto const f4 = ntk.create_node( { d }, 9u );
  auto const f5 = ntk.create_node( { f1 }, 9u );
  auto const f6 = ntk.create_node( { f3, f4, f5 }, 6u );

  ntk.create_po( f6 );
  ntk.create_po( f2 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );
  rinox::analyzers::trackers::arrival_times_tracker<DNtk> tracker( dntk );

  CHECK( tracker.worst_delay() == 3 );

  custom_delay_rewire_params ps;
  rinox::opto::algorithms::delay_resynthesize<Ntk, custom_delay_rewire_params>( ntk, &db, ps );

  CHECK( tracker.worst_delay() == 2 );
}

TEST_CASE( "Delay resynthesis via rewiring - single-output gate with don't cares", "[delay_resynthesis]" )
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
  auto const f1 = ntk.create_node( { a }, 9u );
  auto const f2 = ntk.create_node( { b }, 9u );
  auto const f3 = ntk.create_node( { b, f1 }, 0u );
  auto const f4 = ntk.create_node( { a, f2 }, 0u );
  auto const f5 = ntk.create_node( { f3, f4 }, 1u );
  auto const f6 = ntk.create_node( { a, b }, 1u );
  auto const f7 = ntk.create_node( { f5, f6 }, 0u );
  auto const f8 = ntk.create_node( { a, b }, 8u );

  ntk.create_po( f7 );
  ntk.create_po( f8 );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );
  rinox::analyzers::trackers::arrival_times_tracker<DNtk> tracker( dntk );

  CHECK( tracker.worst_delay() == 4 );

  custom_delay_rewire_params ps;
  ps.window_manager_ps.odc_levels = 3;
  rinox::opto::algorithms::delay_resynthesize<Ntk, custom_delay_rewire_params>( ntk, &db, ps );

  CHECK( tracker.worst_delay() == 2 );
}

TEST_CASE( "Delay resynthesis via rewiring - multiple-output gate without don't cares", "[delay_resynthesis]" )
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
  auto const f1 = ntk.create_node( { c, d, e }, 4u );
  auto const f2 = ntk.create_node( { a, b }, 0u );
  auto const f3 = ntk.create_node( { c, d }, 0u );
  auto const f4 = ntk.create_node( { e, f3 }, 0u );
  auto const f5 = ntk.create_node( { c, f2, f4 }, { 6u, 7u } );

  ntk.create_po( f1 );
  ntk.create_po( { f5.index, 0 } );
  ntk.create_po( { f5.index, 1 } );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );
  rinox::analyzers::trackers::arrival_times_tracker<DNtk> tracker( dntk );

  CHECK( tracker.worst_delay() == 3 );

  custom_delay_rewire_params ps;
  rinox::opto::algorithms::delay_resynthesize<Ntk, custom_delay_rewire_params>( ntk, &db, ps );

  CHECK( tracker.worst_delay() == 2 );
}

TEST_CASE( "Delay resynthesis via rewiring - multiple-output gate with don't cares", "[delay_resynthesis]" )
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
  auto const f1 = ntk.create_node( { a }, 9u );                                // 6
  auto const f2 = ntk.create_node( { b }, 9u );                                // 6
  auto const f3 = ntk.create_node( { a, f2 }, 0u );                            // 6
  auto const f4 = ntk.create_node( { b, f1 }, 0u );                            // 6
  auto const f5 = ntk.create_node( { f3, f4 }, 1u );                           // 6
  auto const f6 = ntk.create_node( { f5, b, c }, { 6u, 7u } );                 // 12
  auto const f7 = ntk.create_node( { a, b }, 8u );                             // 14
  auto const f8 = ntk.create_node( { { f6.index, 0 }, { f6.index, 1 } }, 2u ); // 6
  auto const f9 = ntk.create_node( { f7, f8 }, 0u );                           // 6
  auto const f10 = ntk.create_node( { a, b }, 1u );                            // 6

  ntk.create_po( f9 );
  ntk.create_po( f10 );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );
  rinox::analyzers::trackers::arrival_times_tracker<DNtk> tracker( dntk );

  CHECK( tracker.worst_delay() == 6 );

  custom_delay_rewire_params ps;
  ps.window_manager_ps.odc_levels = 3;
  rinox::opto::algorithms::delay_resynthesize<Ntk, custom_delay_rewire_params>( ntk, &db, ps );

  CHECK( tracker.worst_delay() == 3.0 );
}

struct custom_delay_struct_params1 : rinox::opto::algorithms::default_resynthesis_params<8u>
{
  bool try_rewire = false;
  bool try_struct = true;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  int32_t odc_levels = 0;
  static constexpr uint32_t max_cuts_size = 6u;
  static constexpr uint32_t max_num_vars = 6u;
};

struct custom_delay_struct_params2 : rinox::opto::algorithms::default_resynthesis_params<8u>
{
  bool try_rewire = false;
  bool try_struct = true;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  int32_t odc_levels = 2;
  static constexpr uint32_t max_cuts_size = 6u;
};

TEST_CASE( "Delay resynthesis via cut rewriting - single-output gate without don't cares", "[delay_resynthesis]" )
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
  list.add_output( list.add_gate( { list.add_gate( { 0, 1 }, 8 ), 2 }, 1 ) );
  CHECK( db.add( list ) );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f1 = ntk.create_node( { c, a }, 0u );
  auto const f2 = ntk.create_node( { b }, 9u );
  auto const f3 = ntk.create_node( { f1, f2 }, 0u );
  auto const f4 = ntk.create_node( { f3 }, 9u );

  ntk.create_po( f4 );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );
  rinox::analyzers::trackers::arrival_times_tracker<DNtk> tracker( dntk );

  CHECK( tracker.worst_delay() == 3 );

  custom_delay_struct_params1 ps;
  rinox::opto::algorithms::delay_resynthesize<Ntk, custom_delay_struct_params1>( ntk, &db, ps );

  CHECK( tracker.worst_delay() == 2 );
}

struct custom_delay_struct_params3 : rinox::opto::algorithms::default_resynthesis_params<6u>
{
  bool try_rewire = false;
  bool try_struct = true;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  static constexpr uint32_t max_cuts_size = 3u;
  static constexpr uint32_t max_num_vars = 3u;
};

TEST_CASE( "Delay resynthesis via cut rewriting - single-output gate with don't cares", "[delay_resynthesis]" )
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
  list.add_output( list.add_gate( { list.add_gate( { 0, 1 }, 2 ), 2 }, 0 ) );
  CHECK( db.add( list ) );

  Ntk ntk( gates );
  auto const c = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const a = ntk.create_pi();
  auto const f1 = ntk.create_node( { a, b }, 0u );
  auto const f2 = ntk.create_node( { ntk.create_node( { a }, 9u ), ntk.create_node( { b }, 9u ) }, 8u );
  auto const f3 = ntk.create_node( { f2, c }, 0u );
  auto const f4 = ntk.create_node( { f3, f1 }, 1u );

  ntk.create_po( f4 );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );
  rinox::analyzers::trackers::arrival_times_tracker<DNtk> tracker( dntk );

  CHECK( tracker.worst_delay() == 4 );

  custom_delay_struct_params3 ps;
  ps.window_manager_ps.odc_levels = 3;
  rinox::opto::algorithms::delay_resynthesize<Ntk, custom_delay_struct_params3>( ntk, &db, ps );

  CHECK( tracker.worst_delay() == 3 );
}

struct custom_delay_window_params1 : rinox::opto::algorithms::default_resynthesis_params<6u>
{
  bool try_rewire = false;
  bool try_struct = false;
  bool try_window = true;
  bool try_simula = false;
  bool dynamic_database = false;
  static constexpr uint32_t max_cuts_size = 3u;
  static constexpr uint32_t max_num_vars = 3u;
};

TEST_CASE( "Delay resynthesis via window rewriting - single-output gate without don't cares", "[delay_resynthesis]" )
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
  list.add_output( list.add_gate( { list.add_gate( { 0, 2 }, 2 ), 1 }, 2 ) );
  CHECK( db.add( list ) );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f1 = ntk.create_node( { b, c }, 1u );
  auto const f2 = ntk.create_node( { a, b }, 0u );
  auto const f3 = ntk.create_node( { b, c }, 0u );
  auto const f4 = ntk.create_node( { a, c }, 0u );
  auto const f5 = ntk.create_node( { a, f1 }, 8u );
  auto const f6 = ntk.create_node( { f3 }, 9u );
  auto const f7 = ntk.create_node( { f2, f4 }, 2u );
  auto const f8 = ntk.create_node( { f5, f6 }, 0u );
  auto const f9 = ntk.create_node( { f8 }, 9u );

  ntk.create_po( f3 );
  ntk.create_po( f7 );
  ntk.create_po( f9 );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );

  rinox::analyzers::trackers::arrival_times_tracker<DNtk> tracker( dntk );
  CHECK( tracker.worst_delay() == 4 );

  custom_delay_window_params1 ps;
  rinox::opto::algorithms::delay_resynthesize<Ntk, custom_delay_window_params1>( ntk, &db, ps );

  CHECK( tracker.worst_delay() == 3 );
}

TEST_CASE( "Delay resynthesis via window rewriting - single-output gate with don't cares", "[delay_resynthesis]" )
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
  list.add_output( list.add_gate( { list.add_gate( { 0, 1 }, 2 ), 2 }, 1 ) );
  CHECK( db.add( list ) );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const f1 = ntk.create_node( { a, b }, 0u );
  auto const f2 = ntk.create_node( { c, d }, 0u );
  auto const f3 = ntk.create_node( { b }, 9u );
  auto const f4 = ntk.create_node( { a, f3 }, 8u );
  auto const f5 = ntk.create_node( { f1, b }, 2u );
  auto const f6 = ntk.create_node( { f2, f5 }, 1u );
  auto const f7 = ntk.create_node( { f4, f6 }, 0u );

  ntk.create_po( f2 );
  ntk.create_po( f7 );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );
  rinox::analyzers::trackers::arrival_times_tracker<DNtk> tracker( dntk );
  CHECK( tracker.worst_delay() == 4 );

  custom_delay_window_params1 ps;
  ps.window_manager_ps.odc_levels = 3;
  rinox::opto::algorithms::delay_resynthesize<Ntk, custom_delay_window_params1>( ntk, &db, ps );

  CHECK( tracker.worst_delay() == 3 );
}