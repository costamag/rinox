// SPDX-License-Identifier: MIT
// Tests for rinox boolean/simd_operations.hpp

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <vector>

#include <lorina/genlib.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <rinox/network/network.hpp>

using namespace mockturtle;
using namespace rinox;
using namespace rinox::network;

std::string const test_library = "GATE   inv1    1 O=!a;            PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                 "GATE   inv2    2 O=!a;            PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                 "GATE   nand2   2 O=!(a*b);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   and2    3 O=a*b;           PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
                                 "GATE   xor2    4 O=a^b;           PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                 "GATE   mig3    3 O=a*b+a*c+b*c;   PIN * INV 1 999 2.0 0.2 2.0 0.2\n"
                                 "GATE   xor3    5 O=a^b^c;         PIN * UNKNOWN 2 999 3.0 0.5 3.0 0.5\n"
                                 "GATE   buf     2 O=a;             PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                 "GATE   zero    0 O=CONST0;\n"
                                 "GATE   one     0 O=CONST1;\n"
                                 "GATE   ha      5 C=a*b;           PIN * INV 1 999 1.7 0.4 1.7 0.4\n"
                                 "GATE   ha      5 S=!a*b+a*!b;     PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 C=a*b+a*c+b*c;   PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 S=a^b^c;         PIN * INV 1 999 3.0 0.4 3.0 0.4";

TEST_CASE( "Cell-based Bound network: Primary I / O and constants", "[network]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = bound_network::signal;

  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound_network ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f1 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, 5 );          // maj3
  auto const f2 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { 12, 13 } ); // fa
  signal const carry = signal{ f2.index, 0 };                                             // carry output
  signal const sum = signal{ f2.index, 1 };                                               // sum output
  ntk.create_po( f1 );
  /* use the carry bit as an output */
  ntk.create_po( carry );
  /* create a new node taking the carry signal as input */
  auto const f3 = ntk.create_node<std::vector<signal>, false>( { sum, f1 }, 2u );
  ntk.create_po( f3 );
  CHECK( ntk.is_multioutput( "fa" ) );
  CHECK( !ntk.is_multioutput( "nand2" ) );
  CHECK( ntk.is_combinational() );
  CHECK( ntk.get_constant( true ) == signal{ 1, 0 } );
  CHECK( ntk.get_constant( false ) == signal{ 0, 0 } );
  CHECK( !ntk.is_multioutput( f1.index ) );
  CHECK( ntk.is_multioutput( f2.index ) );
  CHECK( !ntk.is_multioutput( f3.index ) );
  CHECK( !ntk.is_constant( a.index ) );
  CHECK( !ntk.is_constant( f1.index ) );
  CHECK( !ntk.is_constant( f2.index ) );
  CHECK( !ntk.is_constant( f3.index ) );
  CHECK( ( ntk.is_pi( a.index ) && ntk.is_ci( a.index ) ) );
  CHECK( ( ntk.is_pi( b.index ) && ntk.is_ci( b.index ) ) );
  CHECK( ( ntk.is_pi( c.index ) && ntk.is_ci( c.index ) ) );
  CHECK( ( !ntk.is_pi( f1.index ) && !ntk.is_ci( f1.index ) ) );
  CHECK( ( !ntk.is_pi( f2.index ) && !ntk.is_ci( f2.index ) ) );
  CHECK( ( !ntk.is_pi( f3.index ) && !ntk.is_ci( f3.index ) ) );
  CHECK( !ntk.is_po( a ) );
  CHECK( !ntk.is_po( b ) );
  CHECK( !ntk.is_po( c ) );
  CHECK( ntk.is_po( f1 ) );
  // by default, teh signal returned by create node f2 is { f2.index, 0 } = carry
  CHECK( f2 == carry );
  CHECK( ntk.is_po( carry ) );
  CHECK( !ntk.is_po( sum ) );
  CHECK( ntk.is_po( f3 ) );
  CHECK( !ntk.constant_value( 0 ) );
  // any node index different than 0 gives true ( implementation detail )
  CHECK( ntk.constant_value( 1 ) );
  CHECK( ntk.constant_value( 3 ) );
}

TEST_CASE( "Array-based Bound network: Primary I / O and constants", "[network]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 2>;
  using signal = bound_network::signal;

  bound_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  kitty::dynamic_truth_table tt_maj( 3u ), tt_xor( 3u );
  kitty::create_majority( tt_maj );
  kitty::create_parity( tt_xor );
  auto const f1 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, tt_maj );             // maj3
  auto const f2 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { tt_maj, tt_xor } ); // fa
  signal const carry = signal{ f2.index, 0 };                                                     // carry output
  signal const sum = signal{ f2.index, 1 };                                                       // sum output
  ntk.create_po( f1 );
  /* use the carry bit as an output */
  ntk.create_po( carry );
  /* create a new node taking the carry signal as input */
  auto const f3 = ntk.create_node<std::vector<signal>, false>( { sum, f1 }, 2u );
  ntk.create_po( f3 );
  CHECK( ntk.is_combinational() );
  CHECK( ntk.get_constant( true ) == signal{ 1, 0 } );
  CHECK( ntk.get_constant( false ) == signal{ 0, 0 } );
  CHECK( !ntk.is_multioutput( f1.index ) );
  CHECK( ntk.is_multioutput( f2.index ) );
  CHECK( !ntk.is_multioutput( f3.index ) );
  CHECK( !ntk.is_constant( a.index ) );
  CHECK( !ntk.is_constant( f1.index ) );
  CHECK( !ntk.is_constant( f2.index ) );
  CHECK( ( ntk.is_pi( a.index ) && ntk.is_ci( a.index ) ) );
  CHECK( ( ntk.is_pi( b.index ) && ntk.is_ci( b.index ) ) );
  CHECK( ( ntk.is_pi( c.index ) && ntk.is_ci( c.index ) ) );
  CHECK( ( !ntk.is_pi( f1.index ) && !ntk.is_ci( f1.index ) ) );
  CHECK( ( !ntk.is_pi( f2.index ) && !ntk.is_ci( f2.index ) ) );
  CHECK( ( !ntk.is_pi( f3.index ) && !ntk.is_ci( f3.index ) ) );
  CHECK( !ntk.is_po( a ) );
  CHECK( !ntk.is_po( b ) );
  CHECK( !ntk.is_po( c ) );
  CHECK( ntk.is_po( f1 ) );
  // by default, teh signal returned by create node f2 is { f2.index, 0 } = carry
  CHECK( f2 == carry );
  CHECK( ntk.is_po( carry ) );
  CHECK( !ntk.is_po( sum ) );
  CHECK( ntk.is_po( f3 ) );
  CHECK( !ntk.constant_value( 0 ) );
  // any node index different than 0 gives true ( implementation detail )
  CHECK( ntk.constant_value( 1 ) );
  CHECK( ntk.constant_value( 3 ) );
}

TEST_CASE( "Cell-based Bound network: Cloning nodes and networks", "[network]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = bound_network::signal;

  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound_network ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f1 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, 5 );          // maj3
  auto const f2 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { 12, 13 } ); // fa
  signal const carry = signal{ f2.index, 0 };                                             // carry output
  signal const sum = signal{ f2.index, 1 };                                               // sum output
  ntk.create_po( f1 );
  ntk.create_po( carry );
  auto const f3 = ntk.create_node<std::vector<signal>, false>( { sum, f1 }, 2u );
  ntk.create_po( f3 );
  auto ntk2 = ntk.clone();
  CHECK( ntk2.size() == ntk.size() );
  CHECK( ntk2.num_pis() == ntk.num_pis() );
  CHECK( ntk2.num_pos() == ntk.num_pos() );
  CHECK( ntk2.num_gates() == ntk.num_gates() );
  CHECK( ntk2.is_combinational() );

  auto const f4 = ntk2.create_node<std::vector<signal>, false>( { a, b, c }, 6 ); // xor3
  CHECK( ntk2.size() - 1 == ntk.size() );
  CHECK( ntk2.num_gates() - 1 == ntk.num_gates() );

  ntk.clone_node( ntk2, f4.index, { a, b, c } );
  CHECK( ntk2.size() == ntk.size() );
  CHECK( ntk2.num_gates() == ntk.num_gates() );
}

TEST_CASE( "Array-based Bound network: Cloning nodes and networks", "[network]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 2>;
  using signal = bound_network::signal;

  bound_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();

  kitty::dynamic_truth_table tt_maj( 3u ), tt_xor( 3u );
  kitty::create_majority( tt_maj );
  kitty::create_parity( tt_xor );

  auto const f1 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, tt_maj );             // maj3
  auto const f2 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { tt_maj, tt_xor } ); // fa
  signal const carry = signal{ f2.index, 0 };                                                     // carry output
  signal const sum = signal{ f2.index, 1 };                                                       // sum output
  ntk.create_po( f1 );
  ntk.create_po( carry );
  auto const f3 = ntk.create_node<std::vector<signal>, false>( { sum, f1 }, 2u );
  ntk.create_po( f3 );
  auto ntk2 = ntk.clone();
  CHECK( ntk2.size() == ntk.size() );
  CHECK( ntk2.num_pis() == ntk.num_pis() );
  CHECK( ntk2.num_pos() == ntk.num_pos() );
  CHECK( ntk2.num_gates() == ntk.num_gates() );
  CHECK( ntk2.is_combinational() );

  auto const f4 = ntk2.create_node<std::vector<signal>, false>( { a, b, c }, tt_xor ); // xor3
  CHECK( ntk2.size() - 1 == ntk.size() );
  CHECK( ntk2.num_gates() - 1 == ntk.num_gates() );

  ntk.clone_node( ntk2, f4.index, { a, b, c } );
  CHECK( ntk2.size() == ntk.size() );
  CHECK( ntk2.num_gates() == ntk.num_gates() );
}

TEST_CASE( "Cell-based Bound network: Substitute multiple-output node with single-output nodes", "[network]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = bound_network::signal;

  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound_network ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f1 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { 12, 13 } ); // fa
  auto const carry = ntk.create_node<std::vector<signal>, false>( { a, b, c }, 5 );       // maj3
  auto const sum = ntk.create_node<std::vector<signal>, false>( { a, b, c }, 6 );         // xor3
  auto const f2 = ntk.create_node<std::vector<signal>, false>( { signal{ f1.index, 0 },
                                                                 signal{ f1.index, 1 } },
                                                               2u ); // create a new node with carry and sum
  ntk.create_po( signal{ f1.index, 0 } );
  ntk.create_po( f2 );
  ntk.create_po( signal{ f1.index, 1 } );
  ntk.create_po( signal{ f1.index, 0 } );

  CHECK( !ntk.is_po( carry ) );
  CHECK( !ntk.is_po( sum ) );
  CHECK( ntk.is_po( signal{ f1.index, 0 } ) );
  CHECK( ntk.is_po( signal{ f1.index, 1 } ) );
  CHECK( ntk.is_po( f2 ) );
  CHECK( ntk.fanout_size( carry.index ) == 0 );
  CHECK( ntk.fanout_size( sum.index ) == 0 );
  CHECK( ntk.fanout_size( f1.index ) == 5 );
  CHECK( ntk.fanout_size( f2.index ) == 1 );
  CHECK( ntk.size() == 9 );

  ntk.substitute_node( f1.index, std::vector<signal>{ carry, sum } );

  CHECK( ntk.is_po( carry ) );
  CHECK( ntk.is_po( sum ) );
  CHECK( !ntk.is_po( signal{ f1.index, 0 } ) );
  CHECK( !ntk.is_po( signal{ f1.index, 1 } ) );
  CHECK( ntk.is_po( f2 ) );
  CHECK( ntk.fanout_size( carry.index ) == 3 );
  CHECK( ntk.fanout_size( sum.index ) == 2 );
  CHECK( ntk.fanout_size( f1.index ) == 0 );
  CHECK( ntk.fanout_size( f2.index ) == 1 );
  CHECK( ntk.is_dead( f1.index ) );
}

TEST_CASE( "Array-based Bound network: Substitute multiple-output node with single-output nodes", "[network]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 2>;
  using signal = bound_network::signal;

  bound_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();

  kitty::dynamic_truth_table tt_maj( 3u ), tt_xor( 3u );
  kitty::create_majority( tt_maj );
  kitty::create_parity( tt_xor );

  auto const f1 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { tt_maj, tt_xor } ); // fa
  auto const carry = ntk.create_node<std::vector<signal>, false>( { a, b, c }, tt_maj );          // maj3
  auto const sum = ntk.create_node<std::vector<signal>, false>( { a, b, c }, tt_xor );            // xor3
  auto const f2 = ntk.create_node<std::vector<signal>, false>( { signal{ f1.index, 0 },
                                                                 signal{ f1.index, 1 } },
                                                               2u ); // create a new node with carry and sum
  ntk.create_po( signal{ f1.index, 0 } );
  ntk.create_po( f2 );
  ntk.create_po( signal{ f1.index, 1 } );
  ntk.create_po( signal{ f1.index, 0 } );

  CHECK( !ntk.is_po( carry ) );
  CHECK( !ntk.is_po( sum ) );
  CHECK( ntk.is_po( signal{ f1.index, 0 } ) );
  CHECK( ntk.is_po( signal{ f1.index, 1 } ) );
  CHECK( ntk.is_po( f2 ) );
  CHECK( ntk.fanout_size( carry.index ) == 0 );
  CHECK( ntk.fanout_size( sum.index ) == 0 );
  CHECK( ntk.fanout_size( f1.index ) == 5 );
  CHECK( ntk.fanout_size( f2.index ) == 1 );
  CHECK( ntk.size() == 9 );

  ntk.substitute_node( f1.index, std::vector<signal>{ carry, sum } );

  CHECK( ntk.is_po( carry ) );
  CHECK( ntk.is_po( sum ) );
  CHECK( !ntk.is_po( signal{ f1.index, 0 } ) );
  CHECK( !ntk.is_po( signal{ f1.index, 1 } ) );
  CHECK( ntk.is_po( f2 ) );
  CHECK( ntk.fanout_size( carry.index ) == 3 );
  CHECK( ntk.fanout_size( sum.index ) == 2 );
  CHECK( ntk.fanout_size( f1.index ) == 0 );
  CHECK( ntk.fanout_size( f2.index ) == 1 );
  CHECK( ntk.is_dead( f1.index ) );
}

TEST_CASE( "Array-based Bound network: Strashing", "[network]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 2>;
  using signal = typename bound_network::signal;

  bound_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();

  kitty::dynamic_truth_table tt_maj( 3u ), tt_xor( 3u ), tt_nand2( 2u );
  kitty::create_majority( tt_maj );
  kitty::create_parity( tt_xor );
  kitty::create_from_binary_string( tt_nand2, "0111" );

  auto const f1 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { tt_maj, tt_xor } ); // fa
  auto const f2 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { tt_maj, tt_xor } ); // fa
  auto const f3 = ntk.create_node<std::vector<signal>, true>( { a, b, c }, { tt_maj, tt_xor } );  // fa
  auto const f4 = ntk.create_node<std::vector<signal>, false>( { a, b }, tt_nand2 );              // nand2
  auto const f5 = ntk.create_node<std::vector<signal>, false>( { a, b }, tt_nand2 );              // nand2
  auto const f6 = ntk.create_node<std::vector<signal>, true>( { a, b }, tt_nand2 );               // nand2

  CHECK( f2 != f1 );
  CHECK( f3 == f1 );
  CHECK( f3 != f2 );
  CHECK( f5 != f4 );
  CHECK( f6 == f4 );
  CHECK( f6 != f5 );
}

TEST_CASE( "Cell-based Bound network: Strashing", "[network]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename bound_network::signal;

  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound_network ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f1 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { 12, 13 } ); // fa
  auto const f2 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { 12, 13 } ); // fa
  auto const f3 = ntk.create_node<std::vector<signal>, true>( { a, b, c }, { 12, 13 } );  // fa
  auto const f4 = ntk.create_node<std::vector<signal>, false>( { a, b }, 2 );             // nand2
  auto const f5 = ntk.create_node<std::vector<signal>, false>( { a, b }, 2 );             // nand2
  auto const f6 = ntk.create_node<std::vector<signal>, true>( { a, b }, 2 );              // nand2

  CHECK( f2 != f1 );
  CHECK( f3 == f1 );
  CHECK( f3 != f2 );
  CHECK( f5 != f4 );
  CHECK( f6 == f4 );
  CHECK( f6 != f5 );
}

TEST_CASE( "Cell-based Bound network: Simulation", "[network]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename bound_network::signal;

  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound_network ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f1 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { 12, 13 } ); // fa
  auto const f2 = ntk.create_node<std::vector<signal>, false>( { a, b }, 2 );             // nand2

  std::vector<kitty::dynamic_truth_table> tts;
  for ( int i = 0; i < 3; ++i )
  {
    tts.emplace_back( 3u );
    kitty::create_nth_var( tts[i], i );
  }
  tts.push_back( ( tts[0] & tts[1] ) | ( tts[0] & tts[2] ) | ( tts[1] & tts[2] ) ); // maj( a, b, c )
  tts.push_back( ( tts[0] ^ tts[1] ) ^ tts[2] );                                    // a ^ b ^ c
  tts.push_back( ~( tts[0] & tts[1] ) );                                            // nand( a, b )
  std::vector<kitty::dynamic_truth_table const*> sim_ptrs( 3 );
  sim_ptrs[0] = &tts[0];
  sim_ptrs[1] = &tts[1];
  sim_ptrs[2] = &tts[2];
  auto res = ntk.compute( f1.index, sim_ptrs );
  CHECK( res.size() == 2 );
  CHECK( kitty::equal( res[0], tts[3] ) );
  CHECK( kitty::equal( res[1], tts[4] ) );
  sim_ptrs.erase( sim_ptrs.begin() + 2 ); // remove c
  res = ntk.compute( f2.index, sim_ptrs );
  CHECK( res.size() == 1 );
  CHECK( kitty::equal( res[0], tts[5] ) );
}

TEST_CASE( "Array-based Bound network: Simulation", "[network]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 2>;
  using signal = typename bound_network::signal;

  bound_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();

  kitty::dynamic_truth_table tt_maj( 3u ), tt_xor( 3u ), tt_nand2( 2u );
  kitty::create_majority( tt_maj );
  kitty::create_parity( tt_xor );
  kitty::create_from_binary_string( tt_nand2, "0111" );

  auto const f1 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { tt_maj, tt_xor } ); // fa
  auto const f2 = ntk.create_node<std::vector<signal>, false>( { a, b }, tt_nand2 );              // nand2

  std::vector<kitty::dynamic_truth_table> tts;
  for ( int i = 0; i < 3; ++i )
  {
    tts.emplace_back( 3u );
    kitty::create_nth_var( tts[i], i );
  }
  tts.push_back( ( tts[0] & tts[1] ) | ( tts[0] & tts[2] ) | ( tts[1] & tts[2] ) ); // maj( a, b, c )
  tts.push_back( ( tts[0] ^ tts[1] ) ^ tts[2] );                                    // a ^ b ^ c
  tts.push_back( ~( tts[0] & tts[1] ) );                                            // nand( a, b )
  std::vector<kitty::dynamic_truth_table const*> sim_ptrs( 3 );
  sim_ptrs[0] = &tts[0];
  sim_ptrs[1] = &tts[1];
  sim_ptrs[2] = &tts[2];
  auto res = ntk.compute( f1.index, sim_ptrs );
  CHECK( res.size() == 2 );
  CHECK( kitty::equal( res[0], tts[3] ) );
  CHECK( kitty::equal( res[1], tts[4] ) );
  sim_ptrs.erase( sim_ptrs.begin() + 2 ); // remove c
  res = ntk.compute( f2.index, sim_ptrs );
  CHECK( res.size() == 1 );
  CHECK( kitty::equal( res[0], tts[5] ) );
}

TEST_CASE( "Cell-based Bound network wrapped in a depth view", "[network]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = bound_network::signal;

  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  bound_network ntk( gates );
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const f1 = ntk.create_node<std::vector<signal>, false>( { a, b, c }, { 12, 13 } ); // fa
  auto const f2 = ntk.create_node<std::vector<signal>, false>( { signal{ f1.index, 0 },
                                                                 signal{ f1.index, 1 } },
                                                               2u );            // create a new node with carry and sum
  auto const f3 = ntk.create_node<std::vector<signal>, false>( { f2, c }, 2u ); // create a new node with carry and sum
  ntk.create_po( signal{ f1.index, 0 } );
  ntk.create_po( f2 );
  ntk.create_po( signal{ f1.index, 1 } );
  ntk.create_po( signal{ f1.index, 0 } );
  ntk.create_po( f3 );

  depth_view dntk( ntk );
  CHECK( dntk.level( f1.index ) == 1u );
  CHECK( dntk.level( f2.index ) == 2u );
  CHECK( dntk.level( f3.index ) == 3u );
}