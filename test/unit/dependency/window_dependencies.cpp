#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <rinox/dependency/window_dependencies.hpp>
#include <rinox/network/network.hpp>
#include <rinox/windowing/window_manager.hpp>
#include <rinox/windowing/window_simulator.hpp>

std::string const test_library = "GATE   and2    1.0 O=a*b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   or2     1.0 O=a+b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   xor2    1.0 O=a^b;                 PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   or3     1.0 O=a+b+c;               PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   maj3    1.0 O=(a*b)+(b*c)+(a*c);   PIN * INV 1   999 1.0 0.0 1.0 0.0\n"
                                 "GATE   inv1    1.0 O=!(a);                PIN * INV 1   999 1.0 0.0 1.0 0.0";

struct custom_window_params
{
  static constexpr uint32_t max_num_leaves = 6u;
  static constexpr uint32_t max_cuts_size = 6u;
  static constexpr uint32_t max_cube_spfd = 12u;
};

struct window_manager_params : rinox::windowing::default_window_manager_params
{
  static constexpr uint32_t max_num_leaves = 6u;
};

TEST_CASE( "Enumerate dependency cuts", "[window_dependencies]" )
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
  fs.push_back( signal{ 0, 0 } );  // 0
  fs.push_back( signal{ 1, 0 } );  // 1
  fs.push_back( ntk.create_pi() ); // 2
  fs.push_back( ntk.create_pi() ); // 3
  fs.push_back( ntk.create_pi() ); // 4
  fs.push_back( ntk.create_pi() ); // 5

  fs.push_back( ntk.create_node( std::vector<signal>{ fs[2], fs[3] }, 2 ) ); // 6
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[4], fs[5] }, 2 ) ); // 7
  fs.push_back( ntk.create_node( std::vector<signal>{ fs[6], fs[7] }, 2 ) ); // 8

  ntk.create_po( fs[6] );
  ntk.create_po( fs[7] );
  ntk.create_po( fs[8] );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );

  window_manager_params ps;
  ps.odc_levels = 4;
  rinox::windowing::window_manager<DNtk> window( dntk, ps, st );
  CHECK( window.run( dntk.get_node( fs[8] ) ) );
  auto const leaves = window.get_inputs();
  auto const divs = window.get_divisors();
  rinox::windowing::window_simulator<DNtk, custom_window_params::max_num_leaves> sim( dntk );
  sim.run( window );

  rinox::dependency::window_dependencies<DNtk, custom_window_params> dep( dntk );
  dep.run( window, sim );

  std::set<std::set<signal>> sets;
  sets.insert( std::set<signal>( { fs[2], fs[3], fs[4], fs[5] } ) );
  sets.insert( std::set<signal>( { fs[2], fs[3], fs[7] } ) );
  sets.insert( std::set<signal>( { fs[4], fs[5], fs[6] } ) );
  sets.insert( std::set<signal>( { fs[6], fs[7] } ) );
  dep.foreach_cut( [&]( auto const& cut, auto i ) {
    std::set<signal> cut_set;
    for ( auto l : cut.leaves )
      cut_set.insert( l );
    CHECK( sets.find( cut_set ) != sets.end() );
  } );
}

TEST_CASE( "Window dependencies for simple inverter", "[window_dependencies]" )
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
  fs.push_back( signal{ 0, 0 } );  // 0
  fs.push_back( signal{ 1, 0 } );  // 1
  fs.push_back( ntk.create_pi() ); // 2

  fs.push_back( ntk.create_node( std::vector<signal>{ fs[2] }, 5 ) ); // 3

  ntk.create_po( fs[3] );

  using DNtk = mockturtle::depth_view<Ntk>;
  rinox::windowing::window_manager_stats st;
  DNtk dntk( ntk );

  window_manager_params ps;
  ps.odc_levels = 0;
  rinox::windowing::window_manager<DNtk> window( dntk, ps, st );
  CHECK( window.run( dntk.get_node( fs[3] ) ) );
  auto const leaves = window.get_inputs();
  auto const divs = window.get_divisors();
  rinox::windowing::window_simulator<DNtk, custom_window_params::max_num_leaves> sim( dntk );
  sim.run( window );

  rinox::dependency::window_dependencies<DNtk, custom_window_params> dep( dntk );
  dep.run( window, sim );

  std::vector<std::vector<signal>> cuts;
  dep.foreach_cut( [&]( auto const& cut, auto i ) {
    cuts.emplace_back();
    for ( auto l : cut.leaves )
      cuts.back().push_back( l );
  } );
  CHECK( cuts.size() == 1 );
  CHECK( cuts[0].size() == 1 );
  CHECK( cuts[0][0] == fs[2] );
}
