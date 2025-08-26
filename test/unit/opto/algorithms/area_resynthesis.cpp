#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <rinox/network/network.hpp>
#include <rinox/opto/algorithms/resynthesize.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>

std::string const test_library = "GATE   and2    1.0 O=a*b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   or2     1.0 O=a+b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   xor2    0.5 O=a^b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   or3     1.0 O=a+b+c;               PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   and3    1.0 O=((a*b)*c);           PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   maj3    1.0 O=(a*b)+(b*c)+(a*c);   PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   fa      1.0 C=a*b+a*c+b*c;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   fa      1.0 S=a^b^c;               PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   nand2   1.0 O=!(a*b);              PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   inv1    1.0 O=!a;                  PIN * INV 1   999 1.0 0.0 1.0 0.0";

struct custom_area_rewire_params : rinox::opto::algorithms::default_resynthesis_params<8u>
{
  bool try_rewire = true;
  bool try_struct = false;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  int32_t odc_levels = 0;
  static constexpr uint32_t max_cuts_size = 6u;
};

TEST_CASE( "Area resynthesis via rewiring - single-output gate without don't cares", "[area_resynthesis]" )
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
  auto const f5 = ntk.create_node( { f2, f4 }, 0u );

  ntk.create_po( f1 );
  ntk.create_po( f5 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );
  custom_area_rewire_params ps;
  rinox::opto::algorithms::area_resynthesize<DNtk, Db, custom_area_rewire_params>( dntk, db, ps );
  CHECK( ntk.area() == 3 );
}

TEST_CASE( "Area resynthesis via rewiring - single-output gate with don't cares", "[area_resynthesis]" )
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
  auto const f2 = ntk.create_node( { c, d }, 1u );
  auto const f3 = ntk.create_node( { c, d }, 0u );
  auto const f4 = ntk.create_node( { c, d }, 2u );
  auto const f5 = ntk.create_node( { f1, f2 }, 0u );
  auto const f6 = ntk.create_node( { f3, f5 }, 1u );
  auto const f7 = ntk.create_node( { f3, f4 }, 0u );

  ntk.create_po( f6 );
  ntk.create_po( f7 );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );
  custom_area_rewire_params ps;
  ps.window_manager_ps.odc_levels = 3;
  rinox::opto::algorithms::area_resynthesize<DNtk, Db, custom_area_rewire_params>( dntk, db, ps );
  CHECK( ntk.area() == 5.5 );
}

TEST_CASE( "Area resynthesis via rewiring - multiple-output gate without don't cares", "[area_resynthesis]" )
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
  custom_area_rewire_params ps;
  rinox::opto::algorithms::area_resynthesize<DNtk, Db, custom_area_rewire_params>( dntk, db, ps );
  CHECK( ntk.area() == 3 );
}

TEST_CASE( "Area resynthesis via rewiring - multiple-output gate with don't cares", "[area_resynthesis]" )
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
  auto const f1 = ntk.create_node( { a, b }, 0u );
  auto const f2 = ntk.create_node( { c, d }, 1u );
  auto const f3 = ntk.create_node( { c, d }, 0u );
  auto const f4 = ntk.create_node( { c, d }, 2u );
  auto const f5 = ntk.create_node( { e, f1, f2 }, { 6u, 7u } );
  auto const f6 = ntk.create_node( { { f5.index, 0 }, e }, 0u );
  auto const f7 = ntk.create_node( { { f5.index, 1 }, f4 }, 0u );
  auto const f8 = ntk.create_node( { f6, f7 }, 0u );
  auto const f9 = ntk.create_node( { f3, f8 }, 1u );

  ntk.create_po( f9 );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );
  custom_area_rewire_params ps;
  rinox::opto::algorithms::area_resynthesize<DNtk, Db, custom_area_rewire_params>( dntk, db, ps );
  CHECK( ntk.area() == 7.5 );
}


struct custom_area_struct_params1 : rinox::opto::algorithms::default_resynthesis_params<8u>
{
  bool try_rewire = false;
  bool try_struct = true;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  int32_t odc_levels = 0;
  static constexpr uint32_t max_cuts_size = 6u;
};

struct custom_area_struct_params2 : rinox::opto::algorithms::default_resynthesis_params<8u>
{
  bool try_rewire = false;
  bool try_struct = true;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  int32_t odc_levels = 2;
  static constexpr uint32_t max_cuts_size = 6u;
};

TEST_CASE( "Area resynthesis via cut rewriting - single-output gate without don't cares", "[area_resynthesis]" )
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
  custom_area_struct_params1 ps;
  rinox::opto::algorithms::area_resynthesize<DNtk, Db, custom_area_struct_params1>( dntk, db, ps );
  CHECK( ntk.area() == 2.0 );
}

struct custom_area_struct_params3 : rinox::opto::algorithms::default_resynthesis_params<6u>
{
  bool try_rewire = false;
  bool try_struct = true;
  bool try_window = false;
  bool try_simula = false;
  bool dynamic_database = false;
  static constexpr uint32_t max_cuts_size = 3u;
};

TEST_CASE( "Area resynthesis via cut rewriting - single-output gate with don't cares", "[area_resynthesis]" )
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

  rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> list1, list2;
  list1.add_inputs( MaxNumVars );
  list2.add_inputs( MaxNumVars );
  list1.add_output( list1.add_gate( { list1.add_gate( { 0, 1 }, 1 ), 2 }, 0 ) );
  list2.add_output( list2.add_gate( { list2.add_gate( { 0, 1 }, 2 ), 2 }, 0 ) );
  CHECK( db.add( list1 ) );
  CHECK( db.add( list2 ) );

  Ntk ntk( gates );
  auto const c = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const a = ntk.create_pi();
  auto const f1 = ntk.create_node( { a, b }, 0u );
  auto const f2 = ntk.create_node( { a, b }, 1u );
  auto const f3 = ntk.create_node( { f2, c }, 0u );
  auto const f4 = ntk.create_node( { f3, f1 }, 1u );

  ntk.create_po( f4 );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );
  custom_area_struct_params3 ps;
  ps.window_manager_ps.odc_levels = 3;
  rinox::opto::algorithms::area_resynthesize<DNtk, Db, custom_area_struct_params3>( dntk, db, ps );
  CHECK( ntk.area() == 3.5 );
}

struct custom_area_window_params1 : rinox::opto::algorithms::default_resynthesis_params<6u>
{
  bool try_rewire = false;
  bool try_struct = false;
  bool try_window = true;
  bool try_simula = false;
  bool dynamic_database = false;
  static constexpr uint32_t max_cuts_size = 3u;
};

TEST_CASE( "Area resynthesis via window rewriting - single-output gate without don't cares", "[area_resynthesis]" )
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
  list.add_output( list.add_gate( { list.add_gate( { 0, 2 }, 2 ), 1 }, 2 ) );
  CHECK( db.add( list ) );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f1 = ntk.create_node( { a, b }, 0u );
  auto const f2 = ntk.create_node( { a, c }, 0u );
  auto const f3 = ntk.create_node( { b, c }, 0u );
  auto const f4 = ntk.create_node( { f1, f2 }, 1u );
  auto const f5 = ntk.create_node( { f4, f3 }, 1u );

  ntk.create_po( f1 );
  ntk.create_po( f2 );
  ntk.create_po( f3 );
  ntk.create_po( f5 );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );
  custom_area_window_params1 ps;
  rinox::opto::algorithms::area_resynthesize<DNtk, Db, custom_area_window_params1>( dntk, db, ps );
  CHECK( ntk.area() == 4.0 );
}

TEST_CASE( "Area resynthesis via window rewriting - single-output gate with don't cares", "[area_resynthesis]" )
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
  custom_area_window_params1 ps;
  ps.window_manager_ps.odc_levels = 3;
  rinox::opto::algorithms::area_resynthesize<DNtk, Db, custom_area_window_params1>( dntk, db, ps );
  CHECK( ntk.area() == 5.5 );
}

template<uint32_t NumLeaves>
struct ResynParamsCC : rinox::opto::algorithms::default_resynthesis_params<NumLeaves>
{
  static constexpr bool do_strashing = true;
  static constexpr uint32_t num_vars_sign = 12;
  static constexpr uint32_t max_cuts_size = 4u;
  static constexpr uint32_t max_cube_spfd = 10u;

  struct window_manager_params : rinox::windowing::default_window_manager_params
  {
    static constexpr uint32_t max_num_leaves = NumLeaves;
    bool preserve_depth = false;
    int32_t odc_levels = 0u;
    uint32_t skip_fanout_limit_for_divisors = 100u;
    uint32_t max_num_divisors{ 128 };
  };

  window_manager_params window_manager_ps;

};

TEST_CASE( "Delay unconstrained rewiring - corner case 1", "[area_resynthesis]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::ifstream in( std::string(RINOX_SOURCE_DIR) + "/techlib/genlib/asap7.genlib" );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  using Db = rinox::databases::mapped_database<Ntk, 4u>;
  Db db( lib );
  Ntk ntk( gates );

  std::unordered_map<uint32_t, uint64_t> f;
  f[177] = ntk.create_pi();
  f[170] = ntk.create_pi();
  f[176] = ntk.create_pi();
  f[169] = ntk.create_pi();
  f[186] = ntk.create_pi();
  f[188] = ntk.create_pi();
  f[177] = ntk.create_node( { signal{ f[55], 0 } , signal{ f[14], 0 } , signal{ f[26], 0 } , signal{ f[104], 0 } , signal{ f[46], 0 } } , { 34 } ).index;
  ntk.create_po( { f[177], 0 } );
  f[170] = ntk.create_node( { signal{ f[55], 0 } , signal{ f[10], 0 } , signal{ f[24], 0 } , signal{ f[104], 0 } , signal{ f[50], 0 } } , { 34 } ).index;
  ntk.create_po( { f[170], 0 } );
  f[176] = ntk.create_node( { signal{ f[117], 0 } , signal{ f[125], 0 } , signal{ f[154], 0 } , signal{ f[175], 0 } } , { 27 } ).index;
  ntk.create_po( { f[176], 0 } );
  f[169] = ntk.create_node( { signal{ f[117], 0 } , signal{ f[125], 0 } , signal{ f[154], 0 } , signal{ f[122], 0 } } , { 27 } ).index;
  ntk.create_po( { f[169], 0 } );
  f[190] = ntk.create_node( { signal{ f[176], 0 } , signal{ f[177], 0 } } , { 4 } ).index;
  ntk.create_po( { f[190], 0 } );
  f[191] = ntk.create_node( { signal{ f[169], 0 } , signal{ f[170], 0 } } , { 4 } ).index;
  ntk.create_po( { f[191], 0 } );
  f[186] = ntk.create_node( { signal{ f[173], 0 } , signal{ f[172], 0 } } , { 6 } ).index;
  ntk.create_po( { f[186], 0 } );
  f[188] = ntk.create_node( { signal{ f[180], 0 } , signal{ f[178], 0 } } , { 6 } ).index;
  ntk.create_po( { f[188], 0 } );
  f[201] = ntk.create_node( { signal{ f[186], 0 } , signal{ f[190], 0 } , signal{ f[188], 0 } } , { 16 } ).index;
  ntk.create_po( { f[201], 0 } );
  f[200] = ntk.create_node( { signal{ f[191], 0 } } , { 3 } ).index;
  ntk.create_po( { f[200], 0 } );
  f[187] = ntk.create_node( { signal{ f[186], 0 } } , { 3 } ).index;
  f[189] = ntk.create_node( { signal{ f[188], 0 } } , { 3 } ).index;
  f[192] = ntk.create_node( { signal{ f[187], 0 } , signal{ f[189], 0 } , signal{ f[190], 0 } , signal{ f[191], 0 } } , { 21 } ).index;
  ntk.create_po( signal{ f[192], 0 } );
  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk { ntk };
  double area_init = ntk.area(); 
  ResynParamsCC<6u> ps;
  ps.try_rewire = true;
  ps.try_struct = false;
  ps.try_window = false;
  ps.try_simula = false;
  ps.dynamic_database = false;
  ps.fanout_limit = 100;
  ps.stop_after_nth = 7;
  rinox::opto::algorithms::resynthesis_stats st;
  rinox::opto::algorithms::area_resynthesize<DNtk, Db, ResynParamsCC<6u>>( dntk, db, ps, &st );
  double area_end = ntk.area(); 
  CHECK( area_end < area_init );
}