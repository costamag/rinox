#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <rinox/network/network.hpp>
#include <rinox/windowing/window_manager.hpp>

#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>

std::string const test_library = "GATE   inv1    1.0 O=!a ;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   and2    1.0 O=a*b;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   xor2    1.0 O=a^b;         PIN * INV 1   999 3.0 0.0 3.0 0.0";

struct window_manager_params1 : rinox::windowing::default_window_manager_params
{
  static constexpr uint32_t max_num_leaves = 3u;
};

struct window_manager_params2 : rinox::windowing::default_window_manager_params
{
  static constexpr uint32_t max_num_leaves = 8u;
};

TEST_CASE( "Window construction with reconvergent structure", "[window_manager]" )
{

  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  Ntk ntk( gates );

  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  /* Database construction from network */
  signal const a = ntk.create_pi();
  signal const b = ntk.create_pi();
  signal const c = ntk.create_pi();

  fs.push_back( ntk.create_node( { a, b }, 1 ) );           // 0
  fs.push_back( ntk.create_node( { a }, 0 ) );              // 1
  fs.push_back( ntk.create_node( { a, b }, 2 ) );           // 2
  fs.push_back( ntk.create_node( { b, c }, 2 ) );           // 3
  fs.push_back( ntk.create_node( { b, c }, 1 ) );           // 4
  fs.push_back( ntk.create_node( { c }, 0 ) );              // 5
  fs.push_back( ntk.create_node( { b, fs[5] }, 1 ) );       // 6
  fs.push_back( ntk.create_node( { fs[0], fs[1] }, 1 ) );   // 7
  fs.push_back( ntk.create_node( { fs[2], fs[1] }, 2 ) );   // 8
  fs.push_back( ntk.create_node( { fs[2], fs[3] }, 2 ) );   // 9
  fs.push_back( ntk.create_node( { fs[4], fs[3] }, 2 ) );   // 10
  fs.push_back( ntk.create_node( { fs[8], fs[9] }, 2 ) );   // 11
  fs.push_back( ntk.create_node( { fs[9], fs[10] }, 2 ) );  // 12
  fs.push_back( ntk.create_node( { fs[11], fs[12] }, 2 ) ); // 13
  fs.push_back( ntk.create_node( { fs[13], fs[7] }, 2 ) );  // 14
  fs.push_back( ntk.create_node( { fs[13], fs[6] }, 2 ) );  // 15
  fs.push_back( ntk.create_node( { fs[14], fs[1] }, 2 ) );  // 16
  fs.push_back( ntk.create_node( { fs[14], fs[15] }, 2 ) ); // 17
  fs.push_back( ntk.create_node( { fs[5], fs[15] }, 2 ) );  // 18
  fs.push_back( ntk.create_node( { fs[16] }, 0 ) );         // 19
  fs.push_back( ntk.create_node( { fs[17] }, 0 ) );         // 20

  ntk.create_po( fs[19] );
  ntk.create_po( fs[20] );
  ntk.create_po( fs[18] );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );

  window_manager_params1 ps1;
  ps1.odc_levels = 3u;

  rinox::windowing::window_manager<DNtk, window_manager_params1> window( dntk, ps1, st );

  CHECK( window.run( dntk.get_node( fs[13] ) ) );
  CHECK( window.is_valid() );
  auto const mffc = window.get_mffc();
  CHECK( mffc == std::vector<typename Ntk::node_index_t>{ 9, 13, 15, 16, 17, 18 } );
  auto const tfos = window.get_tfos();
  CHECK( tfos == std::vector<typename Ntk::node_index_t>{ 20, 19, 23, 21, 22, 24, 25 } );
  auto const outputs = window.get_outputs();
  CHECK( outputs == std::vector<typename Ntk::signal_t>{ fs[20], fs[19], fs[18] } );
  auto const leaves = window.get_inputs();
  CHECK( leaves == std::vector<typename Ntk::signal_t>( { a, b, c } ) );

  window_manager_params2 ps2;
  ps2.odc_levels = 3u;

  rinox::windowing::window_manager<DNtk, window_manager_params2> window2( dntk, ps2, st );

  CHECK( window2.run( dntk.get_node( fs[13] ) ) );
  CHECK( window2.is_valid() );
  auto const mffc2 = window2.get_mffc();
  CHECK( mffc2 == std::vector<typename Ntk::node_index_t>{ 9, 13, 15, 16, 17, 18 } );
  auto const tfos2 = window2.get_tfos();
  CHECK( tfos2 == std::vector<typename Ntk::node_index_t>{ 20, 19, 23, 21, 22, 24, 25 } );
  auto const outputs2 = window2.get_outputs();
  CHECK( outputs2 == std::vector<typename Ntk::signal_t>( { fs[20], fs[19], fs[18] } ) );
  auto const leaves2 = window2.get_inputs();
  CHECK( leaves2 == std::vector<typename Ntk::signal_t>( { a, b, c } ) );
}

std::string const test_library2 = "GATE   and2    1.0 O=a*b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 0
                                  "GATE   or2     1.0 O=a+b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 1
                                  "GATE   xor2    0.5 O=a^b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 2
                                  "GATE   or3     1.0 O=a+b+c;               PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 3
                                  "GATE   and3    1.0 O=((a*b)*c);           PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 4
                                  "GATE   maj3    1.0 O=(a*b)+(b*c)+(a*c);   PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 5
                                  "GATE   fa      1.0 C=a*b+a*c+b*c;         PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 6
                                  "GATE   fa      1.0 S=a^b^c;               PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 7
                                  "GATE   nand2   1.0 O=!(a*b);              PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 8
                                  "GATE   inv1    1.0 O=!a;                  PIN * INV 1   999 1.0 0.0 1.0 0.0\n" // 9
                                  "GATE   xor3    1.0 O=a^b^c;               PIN * INV 1   999 1.0 0.0 1.0 0.0";  // 10

TEST_CASE( "Reconvergence after multiple-output gate", "[window_manager]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library2 );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

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
  auto const f6 = ntk.create_node( { { f4.index, 0 }, f5 }, 1u );
  auto const f7 = ntk.create_node( { { f4.index, 1 }, f5 }, 1u );
  auto const f8 = ntk.create_node( { a, b }, 1u );
  auto const f9 = ntk.create_node( { c, d }, 2u );
  auto const f10 = ntk.create_node( { f8, f9 }, 0u );

  ntk.create_po( f6 );
  ntk.create_po( f7 );
  ntk.create_po( f10 );

  using DNtk = mockturtle::depth_view<Ntk>;
  DNtk dntk( ntk );
  window_manager_params2 ps;
  rinox::windowing::window_manager_stats st;
  ps.odc_levels = 3u;

  rinox::windowing::window_manager<DNtk, window_manager_params2> window( dntk, ps, st );

  CHECK( window.run( dntk.get_node( f4 ) ) );
  auto const mffc = window.get_mffc();
  auto const tfos = window.get_tfos();
  auto const outputs = window.get_outputs();
  auto const leaves = window.get_inputs();
  CHECK( outputs.size() == 2u );
}

TEST_CASE( "Pathological case 2", "[window_manager]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename Ntk::signal;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library2 );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  Ntk ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
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
  window_manager_params2 ps;
  rinox::windowing::window_manager_stats st;
  ps.odc_levels = 3u;

  rinox::windowing::window_manager<DNtk, window_manager_params2> window( dntk, ps, st );

  CHECK( window.run( dntk.get_node( f1 ) ) );
  auto const mffc = window.get_mffc();
  auto const tfos = window.get_tfos();
  auto const outputs = window.get_outputs();
  auto const leaves = window.get_inputs();
  auto const divs = window.get_divisors();
  CHECK( outputs.size() == 1u );
  CHECK( divs.size() == 4u );
}

std::string const test_library3 = "GATE   and2    1.0 O=a*b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                  "GATE   or2     1.0 O=a+b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                  "GATE   xor2    1.0 O=a^b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                  "GATE   or3     1.0 O=a+b+c;               PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                  "GATE   maj3    1.0 O=(a*b)+(b*c)+(a*c);   PIN * INV 1   999 1.0 0.0 1.0 0.0";

struct window_manager_params3 : rinox::windowing::default_window_manager_params
{
  static constexpr uint32_t max_num_leaves = 16u;
};

TEST_CASE( "Managing exploding number of divisors", "[window_manager]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library3 );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  Ntk ntk( gates );

  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  /* Database construction from network */
  fs.push_back( signal{ 0, 0 } );  // 0
  fs.push_back( signal{ 1, 0 } );  // 1
  fs.push_back( ntk.create_pi() ); // 2
  fs.push_back( ntk.create_pi() ); // 3
  fs.push_back( ntk.create_pi() ); // 4
  fs.push_back( ntk.create_pi() ); // 5
  fs.push_back( ntk.create_pi() ); // 6
  fs.push_back( ntk.create_pi() ); // 7

  fs.push_back( ntk.create_node( std::vector<signal>{ fs[2], fs[3] }, 0 ) );          // 8
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[2], fs[3] }, 1 ) );          // 9
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[2], fs[3] }, 2 ) );          // 10
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[4], fs[5] }, 0 ) );          // 11
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[4], fs[5] }, 1 ) );          // 12
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[4], fs[5] }, 2 ) );          // 13
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[6], fs[7] }, 0 ) );          // 14
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[6], fs[7] }, 1 ) );          // 15
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[6], fs[7] }, 2 ) );          // 16
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[9], fs[12], fs[15] }, 4 ) ); // 17
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[8], fs[17], fs[11] }, 3 ) ); // 18
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[18], fs[14] }, 1 ) );        // 19

  ntk.create_po( fs[19] );
  ntk.create_po( fs[10] );
  ntk.create_po( fs[13] );
  ntk.create_po( fs[16] );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );

  window_manager_params3 ps;
  ps.odc_levels = 4u;

  rinox::windowing::window_manager<DNtk, window_manager_params3> window( dntk, ps, st );
  CHECK( window.run( dntk.get_node( fs[17] ) ) );
  CHECK( ( window.num_divisors() < fs.size() ) );
}

TEST_CASE( "Managing incorrect window construction for simple network", "[window_manager]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library3 );
  auto result = lorina::read_genlib( in, mockturtle::genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  Ntk ntk( gates );

  using signal = typename Ntk::signal;
  std::vector<signal> fs;
  /* Database construction from network */
  fs.push_back( signal{ 0, 0 } );  // 0
  fs.push_back( signal{ 1, 0 } );  // 1
  fs.push_back( ntk.create_pi() ); // 2
  fs.push_back( ntk.create_pi() ); // 3
  fs.push_back( ntk.create_pi() ); // 4

  fs.push_back( ntk.create_node( { fs[2], fs[3] }, 0u ) ); // 5
  fs.push_back( ntk.create_node( { fs[2], fs[4] }, 0u ) ); // 6
  fs.push_back( ntk.create_node( { fs[3], fs[4] }, 0u ) ); // 7
  fs.push_back( ntk.create_node( { fs[5], fs[6] }, 1u ) ); // 8
  fs.push_back( ntk.create_node( { fs[8], fs[7] }, 1u ) ); // 9

  ntk.create_po( fs[5] );
  ntk.create_po( fs[6] );
  ntk.create_po( fs[7] );
  ntk.create_po( fs[9] );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );

  window_manager_params3 ps;
  ps.odc_levels = 0u;

  rinox::windowing::window_manager<DNtk, window_manager_params3> window( dntk, ps, st );
  CHECK( window.run( dntk.get_node( fs[9] ) ) );
  CHECK( window.num_inputs() == 3 );
}
