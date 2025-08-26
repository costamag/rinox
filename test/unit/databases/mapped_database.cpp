#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <sstream>
#include <vector>

#include <lorina/genlib.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <rinox/analyzers/trackers/arrival_times_tracker.hpp>
#include <rinox/databases/mapped_database.hpp>
#include <rinox/evaluation/chains/bound_chain.hpp>
#include <rinox/network/network.hpp>

using namespace mockturtle;

std::string const test_library = "GATE   zero    0 O=CONST0;\n"                                                     // 0
                                 "GATE   one     0 O=CONST1;\n"                                                     // 1
                                 "GATE   inv1    1 O=!a;                      PIN * INV 1 999 0.9 0.3 0.9 0.3\n"    // 2
                                 "GATE   inv2    2 O=!a;                      PIN * INV 2 999 1.0 0.1 1.0 0.1\n"    // 3
                                 "GATE   buf     2 O=a;                       PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n" // 4
                                 "GATE   nand    2 O=!(a*b);                  PIN * INV 1 999 1.0 0.2 1.0 0.2\n"    // 5
                                 "GATE   maj3    8 O=(a*b)+(a*c)+(b*c);       PIN * INV 1 999 3.0 0.4 3.0 0.4\n";   // 6

TEST_CASE( "Adding chains implementing projection to the db", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );
  static constexpr uint32_t MaxNumVars = 4u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );

  rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain0, chain1, chain2, chain3;
  chain0.add_inputs( MaxNumVars );
  auto const a = chain0.pi_at( 0 );
  auto const b = chain0.pi_at( 1 );
  auto const c = chain0.pi_at( 2 );
  auto const d = chain0.pi_at( 3 );
  chain0.add_output( c );

  chain1.add_inputs( MaxNumVars );
  chain1.add_output( a );
  chain2.add_inputs( MaxNumVars );
  chain2.add_output( b );
  chain3.add_inputs( MaxNumVars );
  chain3.add_output( d );

  CHECK( db.add( chain0 ) );
  CHECK( !db.add( chain1 ) );
  CHECK( !db.add( chain2 ) );
  CHECK( !db.add( chain3 ) );
}

std::string const symmetric_library =
    "GATE INV                        1.00  Y=!A;                         \n"
    "    PIN  A  UNKNOWN   1 999    15.00     0.00    15.00     0.00     \n"
    "GATE AND2                       2.00  Y=(A * B);                    \n"
    "    PIN  A  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    20.00     0.00    20.00     0.00     \n"
    "GATE MAJ3                       3.00  Y=(A * B) + (A * C) + (B * C);\n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "GATE ASYM                       3.00  Y=((!A * B) + C);             \n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "GATE AND4                       3.00  Y=((A * B) * (C * D));\n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "    PIN  D  UNKNOWN   1 999    45.00     0.00    25.00     0.00     \n"
    "GATE RND4                       3.00  Y=(((!A * B) + C)^D);         \n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "    PIN  D  UNKNOWN   1 999    65.00     0.00    25.00     0.00     \n"
    "GATE XOR2                       2.00  Y=(A ^ B);                    \n"
    "    PIN  A  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    20.00     0.00    20.00     0.00     \n"
    "GATE FA                       3.00  C=(A * B) + (A * C) + (B * C);  \n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "GATE FA                       3.00  S=( (A ^ B) ^ C );              \n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "GATE RND4_2                     3.00  Y=(((A * B) * C)^D);          \n"
    "    PIN  A  UNKNOWN   1 999    35.00     0.00    35.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    30.00     0.00    30.00     0.00     \n"
    "    PIN  C  UNKNOWN   1 999    25.00     0.00    25.00     0.00     \n"
    "    PIN  D  UNKNOWN   1 999    65.00     0.00    25.00     0.00     \n"
    "GATE OR2                        2.00  Y=(A + B);                    \n"
    "    PIN  A  UNKNOWN   1 999    10.00     0.00    10.00     0.00     \n"
    "    PIN  B  UNKNOWN   1 999    10.00     0.00    10.00     0.00     \n"
    "GATE const0                        0.00  Y=CONST0;                  \n"
    "GATE const1                        0.00  Y=CONST1;                  \n";


TEST_CASE( "Inserting chains with one-input node in mapped databases", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );

  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain;
    chain.add_inputs( MaxNumVars );
    chain.add_output( chain.add_gate( { i }, 0 ) );
    CHECK( !( first ^ db.add( chain ) ) );
    first = false;
  }
}

TEST_CASE( "Inserting chains with two-input node in mapped databases", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );

  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain;
      chain.add_inputs( MaxNumVars );
      chain.add_output( chain.add_gate( { i, j }, 1 ) );
      CHECK( !( first ^ db.add( chain ) ) );
      first = false;
    }
  }
}

TEST_CASE( "Inserting symmetric single-node chains with three inputs in mapped databases", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );
  int cnt = 0;
  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      for ( auto k = 0u; k < MaxNumVars; ++k )
      {
        if ( ( i == k ) || ( j == k ) )
          continue;
        ++cnt;
        if ( cnt > 20 )
          break;
        rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain;
        chain.add_inputs( MaxNumVars );
        chain.add_output( chain.add_gate( { i, j, k }, 2 ) );
        if ( first )
          CHECK( db.add( chain ) );
        else
          CHECK( !db.add( chain ) );
        first = false;
      }
    }
  }
}

TEST_CASE( "Inserting asymmetric single-node chains with three inputs in mapped databases", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );
  int cnt = 0;
  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      for ( auto k = 0u; k < MaxNumVars; ++k )
      {
        if ( ( i == k ) || ( j == k ) )
          continue;
        ++cnt;
        if ( cnt > 20 )
          break;
        rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain;
        chain.add_inputs( MaxNumVars );
        chain.add_output( chain.add_gate( { i, j, k }, 3 ) );
        CHECK( !( first ^ db.add( chain ) ) );
        first = false;
      }
    }
  }
}

TEST_CASE( "Inserting symmetric single-node chains with 4 inputs in mapped databases", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );
  int cnt = 0;
  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      for ( auto k = 0u; k < MaxNumVars; ++k )
      {
        if ( ( i == k ) || ( j == k ) )
          continue;

        for ( auto l = 0u; l < MaxNumVars; ++l )
        {
          if ( ( i == l ) || ( j == l ) || ( k == l ) )
            continue;
          ++cnt;
          if ( cnt > 20 )
            break;
          rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain;
          chain.add_inputs( MaxNumVars );
          chain.add_output( chain.add_gate( { i, j, k, l }, 4 ) );
          CHECK( !( first ^ db.add( chain ) ) );
          first = false;
        }
      }
    }
  }
}

TEST_CASE( "Inserting asymmetric single-node chains with 4 inputs in mapped databases", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );
  int cnt = 0;
  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      for ( auto k = 0u; k < MaxNumVars; ++k )
      {
        if ( ( i == k ) || ( j == k ) )
          continue;

        for ( auto l = 0u; l < MaxNumVars; ++l )
        {
          if ( ( i == l ) || ( j == l ) || ( k == l ) )
            continue;
          ++cnt;
          if ( cnt > 20 )
            break;
          rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain;
          chain.add_inputs( MaxNumVars );
          chain.add_output( chain.add_gate( { i, j, k, l }, 5 ) );
          CHECK( !( first ^ db.add( chain ) ) );
          first = false;
        }
      }
    }
  }
}

TEST_CASE( "Inserting two nodes chain in database", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );

  int cnt = 0;
  bool first = true;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    for ( auto j = 0u; j < MaxNumVars; ++j )
    {
      if ( i == j )
        continue;

      for ( auto k = 0u; k < MaxNumVars; ++k )
      {
        if ( ( i == k ) || ( j == k ) )
          continue;

        for ( auto l = 0u; l < MaxNumVars; ++l )
        {
          if ( ( i == l ) || ( j == l ) || ( k == l ) )
            continue;
          ++cnt;
          if ( cnt > 20 )
            break;
          rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain;
          chain.add_inputs( MaxNumVars );
          auto l0 = chain.add_gate( { i, j }, 1 );
          auto l1 = chain.add_gate( { l0, k, l }, 3 );
          chain.add_output( l1 );
          CHECK( !( first ^ db.add( chain ) ) );
          first = false;
        }
      }
    }
  }
}

TEST_CASE( "Dominant and dominated chains in mapped database", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );
  rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain1, chain2, chain3;
  chain1.add_inputs( MaxNumVars );
  chain2.add_inputs( MaxNumVars );
  chain3.add_inputs( MaxNumVars );
  auto const l1_1 = chain1.add_gate( { 1 }, 0 );
  auto const l1_2 = chain1.add_gate( { 5 }, 0 );
  auto const l1_3 = chain1.add_gate( { l1_1, 5 }, 1 );
  auto const l1_4 = chain1.add_gate( { l1_2, 1 }, 1 );
  auto const l1_5 = chain1.add_gate( { l1_3, l1_4 }, 6 );
  chain1.add_output( l1_5 );

  auto const l2_1 = chain2.add_gate( { 4, 0 }, 6 );
  chain2.add_output( l2_1 );

  chain3 = chain1;

  CHECK( db.size() == 0 );
  CHECK( db.num_rows() == 0 );
  CHECK( db.add( chain1 ) );
  CHECK( db.num_rows() == 1 );
  CHECK( db.size() == 1 );
  CHECK( db.add( chain2 ) );
  CHECK( db.size() == 1 );
  CHECK( db.num_rows() == 1 );
  CHECK( !db.add( chain3 ) );
  CHECK( db.num_rows() == 1 );
  CHECK( db.size() == 1 );
}

TEST_CASE( "Saving a mapped database", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 6u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );
  rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain1, chain2, chain3, chain4;
  chain1.add_inputs( MaxNumVars );
  chain2.add_inputs( MaxNumVars );
  chain3.add_inputs( MaxNumVars );
  chain4.add_inputs( MaxNumVars );
  auto const l1_1 = chain1.add_gate( { 1 }, 0 );
  auto const l1_2 = chain1.add_gate( { 5 }, 0 );
  auto const l1_3 = chain1.add_gate( { l1_1, 5 }, 1 );
  auto const l1_4 = chain1.add_gate( { l1_2, 1 }, 1 );
  auto const l1_5 = chain1.add_gate( { l1_3, l1_4 }, 6 );
  chain1.add_output( l1_5 );

  auto const l2_1 = chain2.add_gate( { 4, 0 }, 6 );
  chain2.add_output( l2_1 );

  auto const l3_1 = chain3.add_gate( { 1, 5, 2, 0 }, 4 );
  auto const l3_2 = chain3.add_gate( { l3_1 }, 0 );
  auto const l3_3 = chain3.add_gate( { 3, l3_2 }, 1 );
  chain3.add_output( l3_3 );

  auto const l4_1 = chain4.add_gate( { 2, 0, 3, 1 }, 4 );
  auto const l4_2 = chain4.add_gate( { l4_1 }, 0 );
  auto const l4_3 = chain4.add_gate( { 4, l4_2 }, 1 );
  chain4.add_output( l4_3 );

  CHECK( db.add( chain1 ) );
  CHECK( db.num_rows() == 1 );
  CHECK( db.size() == 1 );
  CHECK( db.add( chain2 ) );
  CHECK( db.num_rows() == 1 );
  CHECK( db.size() == 1 );
  CHECK( db.add( chain3 ) );
  CHECK( db.num_rows() == 2 );
  CHECK( !db.add( chain4 ) );
  CHECK( db.num_rows() == 2 );
  CHECK( db.size() == 2 );
  std::stringstream out;
  db.commit( out );
  std::string const expected =
"module top( x0 , x1 , x2 , x3 , x4 , x5 , y0 , y1 );\n"
"  input x0 ;\n"
"  input x1 ;\n"
"  input x2 ;\n"
"  input x3 ;\n"
"  input x4 ;\n"
"  input x5 ;\n"
"  output y0 ;\n"
"  output y1 ;\n"
"  wire n0 , n1 , n14 , n15 ;\n"
"  const0 g00( .Y (n0) );\n"
"  const1 g01( .Y (n1) );\n"
"  XOR2   g02( .A (x4), .B (x5), .Y (y0) );\n"
"  AND4   g03( .A (x3), .B (x4), .C (x5), .D (x2), .Y (n14) );\n"
"  INV    g04( .A (n14), .Y (n15) );\n"
"  AND2   g05( .A (x0), .B (n15), .Y (y1) );\n"
"endmodule\n";
  CHECK( out.str() == expected );
}

TEST_CASE( "Database look-up with 2-input completely specified function", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename bound_network::signal;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 3u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );

  rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain;
  chain.add_inputs( MaxNumVars );
  chain.add_output( chain.add_gate( { chain.add_gate( { 0 }, 0 ), 1 }, 1 ) );
  CHECK( db.add( chain ) );

  using TT = kitty::static_truth_table<MaxNumVars>;

  bound_network ntk( gates );
  std::array<TT, 3u> xs;
  std::vector<signal> fs;
  std::vector<TT const*> sim_ptrs;
  for ( auto i = 0u; i < 3u; ++i )
  {
    kitty::create_nth_var( xs[i], i );
    sim_ptrs.push_back( &xs[i] );
    fs.push_back( ntk.create_pi() );
  }

  std::vector<uint32_t> input{ 0, 1 };
  std::vector<double> times{ 0, 10 };

  // symmetric function

  rinox::evaluation::chain_simulator<rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED>, TT> sim( lib );

  auto tt = ( ~xs[0] ) & xs[1];
  auto fs_c = fs;
  auto times_c = times;
  auto match = db.boolean_matching( tt, times_c, fs_c );
  CHECK( match );
  db.foreach_entry( *match, [&]( auto const& entry ) {
    auto const n = db.write( entry, ntk, fs_c );

    rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain_res( MaxNumVars );
    rinox::evaluation::chains::extract( chain_res, ntk, fs, ntk.make_signal( n ) );
    sim( chain_res, sim_ptrs );

    auto const res = sim.get_simulation( chain_res, sim_ptrs, chain_res.po_at( 0 ) );
    if ( !kitty::equal( res, tt ) )
    {
      kitty::print_binary( tt );
      std::cout << std::endl;
      kitty::print_binary( res );
      std::cout << std::endl;
    }
    CHECK( kitty::equal( res, tt ) );
  } );
}

TEST_CASE( "Database look-up with 3-input completely specified function", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using signal = typename bound_network::signal;
  std::vector<gate> gates;

  std::istringstream in( symmetric_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );

  static constexpr uint32_t MaxNumVars = 4u;
  rinox::databases::mapped_database<bound_network, MaxNumVars> db( lib );

  rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain;
  chain.add_inputs( MaxNumVars );
  auto const l0 = chain.add_gate( { 0u }, 0u );
  auto const l1 = chain.add_gate( { 1u, 2u }, 1u );
  auto const l2 = chain.add_gate( { l0, l1 }, 10u );
  chain.add_output( l2 );
  CHECK( db.add( chain ) );

  using TT = kitty::static_truth_table<MaxNumVars>;

  bound_network ntk( gates );
  std::array<TT, MaxNumVars> xs;
  std::vector<signal> fs;
  std::vector<TT const*> sim_ptrs;
  for ( auto i = 0u; i < MaxNumVars; ++i )
  {
    kitty::create_nth_var( xs[i], i );
    sim_ptrs.push_back( &xs[i] );
    fs.push_back( ntk.create_pi() );
  }

  std::vector<double> times{ 0, 10, 20, 40 };
  rinox::analyzers::trackers::arrival_times_tracker arrival( ntk, times );

  rinox::evaluation::chain_simulator<rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED>, TT> sim( lib );

  int i = 0;
  auto tt = ( ~xs[1] ) | ( xs[2] & xs[3] );
  auto fs_c = fs;
  auto match = db.boolean_matching( tt, times, fs_c );
  CHECK( match );
  db.foreach_entry( *match, [&]( auto const& entry ) {
    auto const n = db.write( entry, ntk, fs_c );
    ntk.create_po( ntk.make_signal( n ) );
    rinox::evaluation::chains::bound_chain<rinox::network::design_type_t::CELL_BASED> chain_res( MaxNumVars );
    rinox::evaluation::chains::extract( chain_res, ntk, fs, ntk.make_signal( n ) );
    sim( chain_res, sim_ptrs );

    auto const res = sim.get_simulation( chain_res, sim_ptrs, chain_res.po_at( 0 ) );
    CHECK( kitty::equal( res, tt ) );
    CHECK( fs_c[0] == fs[2] );
    CHECK( fs_c[1] == fs[3] );
    CHECK( fs_c[2] == fs[0] );
    CHECK( times[0] == 20 );
    CHECK( times[1] == 40 );
  } );

  CHECK( arrival.worst_delay() == 70 );
}

TEST_CASE( "Loading a mapped database", "[mapped_database]" )
{
  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<gate> gates;

  std::istringstream in_lib( symmetric_library );
  auto result = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( gates );
  static constexpr uint32_t MaxNumVars = 4u;

  std::string const stored =
      "module top( x0 , x1 , x2 , x3 , y0 , y1 );\n"
      "  input x0 ;\n"
      "  input x1 ;\n"
      "  input x2 ;\n"
      "  input x3 ;\n"
      "  output y0 ;\n"
      "  output y1 ;\n"
      "  wire n0 , n1 , n7 , n8 ;\n"
      "  const0 g0( .Y (n0) );\n"
      "  const1 g1( .Y (n1) );\n"
      "  XOR2   g2( .A (x2), .B (x3), .Y (y0) );\n"
      "  INV    g3( .A (x3), .Y (n7) );\n"
      "  AND2   g4( .A (x2), .B (n7), .Y (n8) );\n"
      "  MAJ3   g5( .A (x0), .B (x1), .C (n8), .Y (y1) );\n"
      "endmodule\n";

  rinox::databases::mapped_database<bound_network, 4u> db( lib );
  std::istringstream in( stored );
  db.load( in );

  std::ostringstream out;
  db.commit( out );
  CHECK( out.str() == stored );
}