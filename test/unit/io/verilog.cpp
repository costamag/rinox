#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <string>

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/networks/buffered.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/muxig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/xmg.hpp>
#include <rinox/io/utils/reader.hpp>
#include <rinox/io/verilog/verilog.hpp>
#include <rinox/network/network.hpp>

#include <kitty/kitty.hpp>
#include <lorina/verilog.hpp>

using namespace mockturtle;

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

TEST_CASE( "Read structural verilog to mapped network", "[verilog_reader]" )
{

  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in_lib( test_library );
  auto result_lib = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result_lib == lorina::return_code::success );

  std::string file{
      "module top( x0 , x1 , x2 , y0 , y1 , y2, y3 );\n"
      "  input x0 , x1, x2 ;\n"
      "  output y0 , y1 , y2, y3 ;\n"
      "  wire n4 , n5 , n6 ;\n"
      "  inv1 g0( .a (x0), .O (n4) );\n"
      "  fa   g1( .a (n4), .b (x1), .c (x2), .C (n5), .S (n6) );\n"
      "  inv1 g2( .a (n4), .O (y0) );\n"
      "  xor2 g3( .a (n6), .b (x2), .O (y1) );\n"
      "  buf g4( .a (n5), .O (y2) );\n"
      "  buf g5( .a (n6), .O (y3) );\n"
      "endmodule\n" };

  std::istringstream in_ntk( file );

  bound_network ntk( gates );
  const auto result_ntk = rinox::io::verilog::read_verilog( in_ntk, rinox::io::reader( ntk ) );

  /* structural checks */
  CHECK( result_ntk == lorina::return_code::success );
  CHECK( ntk.num_pis() == 3 );
  CHECK( ntk.num_pos() == 4 );
  CHECK( ntk.size() == 11 ); // 2 constants, 2 PIs, 3 buffers, 1 inverter, 3 crossings, 3 gates

  CHECK( ntk.num_gates() == 6 );

  std::ostringstream out;
  rinox::io::verilog::write_verilog( ntk, out );

  std::string expected = "module top( x0 , x1 , x2 , y0 , y1 , y2 , y3 );\n"
                         "  input x0 ;\n"
                         "  input x1 ;\n"
                         "  input x2 ;\n"
                         "  output y0 ;\n"
                         "  output y1 ;\n"
                         "  output y2 ;\n"
                         "  output y3 ;\n"
                         "  wire n0 , n1 , n5 , n6_0 , n6_1 ;\n"
                         "  zero  g00( .O (n0) );\n"
                         "  one   g01( .O (n1) );\n"
                         "  inv1  g02( .a (x0), .O (n5) );\n"
                         "  inv1  g03( .a (n5), .O (y0) );\n"
                         "  fa    g04( .a (n5), .b (x1), .c (x2), .C (n6_0), .S (n6_1) );\n"
                         "  xor2  g05( .a (n6_1), .b (x2), .O (y1) );\n"
                         "  buf   g06( .a (n6_0), .O (y2) );\n"
                         "  buf   g07( .a (n6_1), .O (y3) );\n"
                         "endmodule\n";
  CHECK( out.str() == expected );
}

TEST_CASE( "Read structural verilog to mapped network  with the inputs permutated", "[verilog_reader]" )
{

  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in_lib( test_library );
  auto result_lib = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result_lib == lorina::return_code::success );

  std::string file{
      "module top( x0 , x1 , x2 , y0 , y1 , y2, y3 );\n"
      "  input x0 , x1, x2 ;\n"
      "  output y0 , y1 , y2, y3 ;\n"
      "  wire n4 , n5 , n6 ;\n"
      "  inv1 g0( .a (x0), .O (n4) );\n"
      "  fa   g1( .b (x1), .C (n5), .a (n4), .c (x2), .S (n6) );\n"
      "  inv1 g2( .a (n4), .O (y0) );\n"
      "  xor2 g3( .a (n6), .O (y1), .b (x2) );\n"
      "  buf g4( .a (n5), .O (y2) );\n"
      "  buf g5( .O (y3), .a (n6) );\n"
      "endmodule\n" };

  std::istringstream in_ntk( file );

  bound_network ntk( gates );
  const auto result_ntk = rinox::io::verilog::read_verilog( in_ntk, rinox::io::reader( ntk ) );

  /* structural checks */
  CHECK( result_ntk == lorina::return_code::success );
  CHECK( ntk.num_pis() == 3 );
  CHECK( ntk.num_pos() == 4 );
  CHECK( ntk.size() == 11 ); // 2 constants, 2 PIs, 3 buffers, 1 inverter, 3 crossings, 3 gates

  CHECK( ntk.num_gates() == 6 );

  std::ostringstream out;
  rinox::io::verilog::write_verilog( ntk, out );

  std::string expected = "module top( x0 , x1 , x2 , y0 , y1 , y2 , y3 );\n"
                         "  input x0 ;\n"
                         "  input x1 ;\n"
                         "  input x2 ;\n"
                         "  output y0 ;\n"
                         "  output y1 ;\n"
                         "  output y2 ;\n"
                         "  output y3 ;\n"
                         "  wire n0 , n1 , n5 , n6_0 , n6_1 ;\n"
                         "  zero  g00( .O (n0) );\n"
                         "  one   g01( .O (n1) );\n"
                         "  inv1  g02( .a (x0), .O (n5) );\n"
                         "  inv1  g03( .a (n5), .O (y0) );\n"
                         "  fa    g04( .a (n5), .b (x1), .c (x2), .C (n6_0), .S (n6_1) );\n"
                         "  xor2  g05( .a (n6_1), .b (x2), .O (y1) );\n"
                         "  buf   g06( .a (n6_0), .O (y2) );\n"
                         "  buf   g07( .a (n6_1), .O (y3) );\n"
                         "endmodule\n";
  CHECK( out.str() == expected );
}

TEST_CASE( "Read structural verilog to mapped network in the format used by abc", "[verilog_reader]" )
{

  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in_lib( test_library );
  auto result_lib = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result_lib == lorina::return_code::success );

  std::string file = R"(module top( \a[0] , \a[1] , \b[0] , \f[0] , \f[1] , \f[2], \f[3] );
    input \a[0] , \a[1] , \b[0] ;
    output \f[0] , \f[1] , \f[2], \f[3] ;
    wire n4 , n5 , n6 ;
    inv1 g0( .a (\a[0]), .O (n4) );
    fa   g1( .b (\a[1]), .S (n5), .a (n4), .c (\b[0]), .C (n6) );
    inv1 g2( .a (n4), .O (\f[0]) );
    xor2 g3( .a (n6), .O (\f[1]), .b (\b[0]) );
    buf g4( .a (n5), .O (\f[2]) );
    buf g5( .O (\f[3]), .a (n6) );
    endmodule)";

  std::istringstream in_ntk( file );

  bound_network ntk( gates );
  const auto result_ntk = rinox::io::verilog::read_verilog( in_ntk, rinox::io::reader( ntk ) );

  /* structural checks */
  CHECK( result_ntk == lorina::return_code::success );
  CHECK( ntk.num_pis() == 3 );
  CHECK( ntk.num_pos() == 4 );
  CHECK( ntk.size() == 11 ); // 2 constants, 2 PIs, 3 buffers, 1 inverter, 3 crossings, 3 gates

  CHECK( ntk.num_gates() == 6 );

  std::ostringstream out;
  rinox::io::verilog::write_verilog( ntk, out );

  std::string expected =
      "module top( a , b , f );\n"
      "  input [0:1] a ;\n"
      "  input b ;\n"
      "  output [0:3] f ;\n"
      "  wire n0 , n1 , n5 , n6_0 , n6_1 ;\n"
      "  zero  g00( .O (n0) );\n"
      "  one   g01( .O (n1) );\n"
      "  inv1  g02( .a (a[0]), .O (n5) );\n"
      "  inv1  g03( .a (n5), .O (f[0]) );\n"
      "  fa    g04( .a (n5), .b (a[1]), .c (b), .S (n6_0), .C (n6_1) );\n"
      "  xor2  g05( .a (n6_1), .b (b), .O (f[1]) );\n"
      "  buf   g06( .a (n6_0), .O (f[2]) );\n"
      "  buf   g07( .a (n6_1), .O (f[3]) );\n"
      "endmodule\n";

  CHECK( out.str() == expected );
}

TEST_CASE( "Corner case 1 - assign constant", "[verilog_reader]" )
{

  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in_lib( test_library );
  auto result_lib = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result_lib == lorina::return_code::success );

  std::string file =
      R"(module top(a, y);
  input a;
  output y;
  wire _7_;

  assign _7_ = 1'h0;

  buf g0 (.a(_7_), .O(y));
  endmodule
  )";

  std::istringstream in_ntk( file );

  bound_network ntk( gates );
  const auto result_ntk = rinox::io::verilog::read_verilog( in_ntk, rinox::io::reader( ntk ) );

  /* structural checks */
  CHECK( result_ntk == lorina::return_code::success );
  CHECK( ntk.num_pis() == 1 );
  CHECK( ntk.num_pos() == 1 );
  CHECK( ntk.size() == 4 ); // 2 constants, 2 PIs, 3 buffers, 1 inverter, 3 crossings, 3 gates

  CHECK( ntk.num_gates() == 1 );

  std::ostringstream out;
  rinox::io::verilog::write_verilog( ntk, out );

  std::string expected =
      "module top( a , y );\n"
      "  input a ;\n"
      "  output y ;\n"
      "  wire n0 , n1 ;\n"
      "  zero  g0( .O (n0) );\n"
      "  one   g1( .O (n1) );\n"
      "  buf   g2( .a (n0), .O (y) );\n"
      "endmodule\n";

  CHECK( out.str() == expected );
}

TEST_CASE( "Corner case 2 - vector assignment with positive polarity", "[verilog_reader]" )
{

  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in_lib( test_library );
  auto result_lib = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result_lib == lorina::return_code::success );

  std::string file =
      R"(module top( \xyz_inst.a[0] , \xyz_inst.a[1] , \xyz_inst.a[2] , \xyz_inst.a[3] ,
            \xyz_inst.i0.a[0] , \xyz_inst.i0.a[1] , \xyz_inst.i0.a[2] , \xyz_inst.i0.a[3] );
    // Port directions
    input  \xyz_inst.a[0] , \xyz_inst.a[1] , \xyz_inst.a[2] , \xyz_inst.a[3] ;
    output \xyz_inst.i0.a[0] , \xyz_inst.i0.a[1] , \xyz_inst.i0.a[2] , \xyz_inst.i0.a[3] ;

    // Intermediate net for variety
    wire tmp;

    // Vector-to-vector assign with escaped IDs
    assign { \xyz_inst.i0.a[3] , \xyz_inst.i0.a[2] , \xyz_inst.i0.a[1] , tmp } =
           { \xyz_inst.a[3]    , \xyz_inst.a[2]    , \xyz_inst.a[1]    , \xyz_inst.a[0]    };

    // Some trivial logic to touch tmp
    buf g1 (.a(tmp), .O(\xyz_inst.i0.a[0]));

endmodule
  )";

  std::istringstream in_ntk( file );

  bound_network ntk( gates );
  const auto result_ntk = rinox::io::verilog::read_verilog( in_ntk, rinox::io::reader( ntk ) );

  /* structural checks */
  CHECK( result_ntk == lorina::return_code::success );
  CHECK( ntk.num_pis() == 4 );
  CHECK( ntk.num_pos() == 4 );
  CHECK( ntk.size() == 7 ); // 2 constants, 4 PIs, 2 buffers
  CHECK( ntk.num_gates() == 1 );

  std::ostringstream out;
  rinox::io::verilog::write_verilog( ntk, out );

  std::string expected =
      "module top( xyz_inst.a , xyz_inst.i0.a );\n"
      "  input [0:3] xyz_inst.a ;\n"
      "  output [0:3] xyz_inst.i0.a ;\n"
      "  wire n0 , n1 ;\n"
      "  zero  g0( .O (n0) );\n"
      "  one   g1( .O (n1) );\n"
      "  buf   g2( .a (xyz_inst.a[0]), .O (xyz_inst.i0.a[0]) );\n"
      "  assign xyz_inst.i0.a[1] = xyz_inst.a[1] ;\n"
      "  assign xyz_inst.i0.a[2] = xyz_inst.a[2] ;\n"
      "  assign xyz_inst.i0.a[3] = xyz_inst.a[3] ;\n"
      "endmodule\n";

  CHECK( out.str() == expected );
}

TEST_CASE( "Corner case 3 - vector assignment with explicit positive polarity", "[verilog_reader]" )
{

  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in_lib( test_library );
  auto result_lib = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result_lib == lorina::return_code::success );

  std::string file =
      R"(module top( \xyz_inst.f00[0] , \xyz_inst.f01[0] , \xyz_inst.f02[0] , \xyz_inst.f03[0] ,
            \xyz_inst.f04[0] , \xyz_inst.f05[0] , \xyz_inst.f06[0] , \xyz_inst.f07[0] ,
            \xyz_inst.out0[0], \xyz_inst.out0[1], \xyz_inst.out0[2], \xyz_inst.out0[3],
            \xyz_inst.out0[4], \xyz_inst.out0[5], \xyz_inst.out0[6], \xyz_inst.out0[7] );

    // Port directions
    input  \xyz_inst.f00[0] , \xyz_inst.f01[0] , \xyz_inst.f02[0] , \xyz_inst.f03[0] ,
           \xyz_inst.f04[0] , \xyz_inst.f05[0] , \xyz_inst.f06[0] , \xyz_inst.f07[0] ;
    output \xyz_inst.out0[0], \xyz_inst.out0[1], \xyz_inst.out0[2], \xyz_inst.out0[3],
           \xyz_inst.out0[4], \xyz_inst.out0[5], \xyz_inst.out0[6], \xyz_inst.out0[7] ;

    // Internal nets for a couple of small gates
    wire n1, n2;

    // Two trivial logic gates
    and2 g0 (.a(\xyz_inst.f00[0]), .b(\xyz_inst.f01[0]), .O(n1));
    xor2  g1 (.a(\xyz_inst.f02[0]), .b(n1),               .O(n2));

    // Weird continuous assignment: concatenations with an addition in between
    assign { \xyz_inst.out0[7], \xyz_inst.out0[6], \xyz_inst.out0[5], \xyz_inst.out0[4],
             \xyz_inst.out0[3], \xyz_inst.out0[2], \xyz_inst.out0[1], \xyz_inst.out0[0] }
         = + { \xyz_inst.f07[0], \xyz_inst.f06[0], \xyz_inst.f05[0], \xyz_inst.f04[0],
                \xyz_inst.f03[0], n2, \xyz_inst.f01[0], \xyz_inst.f00[0] };

endmodule)";

  std::istringstream in_ntk( file );

  bound_network ntk( gates );
  const auto result_ntk = rinox::io::verilog::read_verilog( in_ntk, rinox::io::reader( ntk ) );

  /* structural checks */
  CHECK( result_ntk == lorina::return_code::success );
  CHECK( ntk.num_pis() == 8 );
  CHECK( ntk.num_pos() == 8 );
  CHECK( ntk.size() == 12 ); // 2 constants, 4 PIs, 2 buffers
  CHECK( ntk.num_gates() == 2 );

  std::ostringstream out;
  rinox::io::verilog::write_verilog( ntk, out );

  std::string expected =
      "module top( xyz_inst.f00 , xyz_inst.f01 , xyz_inst.f02 , xyz_inst.f03 , xyz_inst.f04 , xyz_inst.f05 , xyz_inst.f06 , xyz_inst.f07 , xyz_inst.out0 );\n"
      "  input xyz_inst.f00 ;\n"
      "  input xyz_inst.f01 ;\n"
      "  input xyz_inst.f02 ;\n"
      "  input xyz_inst.f03 ;\n"
      "  input xyz_inst.f04 ;\n"
      "  input xyz_inst.f05 ;\n"
      "  input xyz_inst.f06 ;\n"
      "  input xyz_inst.f07 ;\n"
      "  output [0:7] xyz_inst.out0 ;\n"
      "  wire n0 , n1 , n10 ;\n"
      "  zero  g00( .O (n0) );\n"
      "  one   g01( .O (n1) );\n"
      "  and2  g02( .a (xyz_inst.f00), .b (xyz_inst.f01), .O (n10) );\n"
      "  xor2  g03( .a (xyz_inst.f02), .b (n10), .O (xyz_inst.out0[2]) );\n"
      "  assign xyz_inst.out0[0] = xyz_inst.f00 ;\n"
      "  assign xyz_inst.out0[1] = xyz_inst.f01 ;\n"
      "  assign xyz_inst.out0[3] = xyz_inst.f03 ;\n"
      "  assign xyz_inst.out0[4] = xyz_inst.f04 ;\n"
      "  assign xyz_inst.out0[5] = xyz_inst.f05 ;\n"
      "  assign xyz_inst.out0[6] = xyz_inst.f06 ;\n"
      "  assign xyz_inst.out0[7] = xyz_inst.f07 ;\n"
      "endmodule\n";

  CHECK( out.str() == expected );
}

TEST_CASE( "Corner case 4 - vector assignment with explicit negative polarity", "[verilog_reader]" )
{

  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in_lib( test_library );
  auto result_lib = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result_lib == lorina::return_code::success );

  std::string file =
      R"(module top( \a[0] , \a[1] , \b[0] , \b[1] ,
            \xyz_inst.out1[0] , \xyz_inst.out1[1] , \xyz_inst.out1[2] , \xyz_inst.out1[3] ,
            \xyz_inst.out1[4] , \xyz_inst.out1[5] , \xyz_inst.out1[6] , \xyz_inst.out1[7] );

    // Ports
    input  \a[0] , \a[1] , \b[0] , \b[1] ;
    output \xyz_inst.out1[0] , \xyz_inst.out1[1] , \xyz_inst.out1[2] , \xyz_inst.out1[3] ,
           \xyz_inst.out1[4] , \xyz_inst.out1[5] , \xyz_inst.out1[6] , \xyz_inst.out1[7] ;

    // Internal nets (all are used)
    wire \xyz_inst.f07[1] , \xyz_inst.f06[1] , \xyz_inst.f05[1] , \xyz_inst.f04[1] , \xyz_inst.f03[1] ;
    wire t0;

    // A couple of gates to drive the f** nets
    and2 g0 (.a(\a[0]), .b(\b[0]), .O(\xyz_inst.f07[1]));
    xor2 g1 (.a(\a[1]), .b(\b[1]), .O(\xyz_inst.f06[1]));
    inv1 g2 (.a(\a[0]), .O(t0));
    xor2 g3 (.a(t0), .b(\b[1]), .O(\xyz_inst.f05[1]));
    buf  g4 (.a(\a[1]), .O(\xyz_inst.f04[1]));
    buf  g5 (.a(\b[0]), .O(\xyz_inst.f03[1]));

    // Subtraction with concatenations (8-bit wide on each side)
    assign { \xyz_inst.out1[7] , \xyz_inst.out1[6] , \xyz_inst.out1[5] , \xyz_inst.out1[4] ,
             \xyz_inst.out1[3] , \xyz_inst.out1[2] , \xyz_inst.out1[1] , \xyz_inst.out1[0] } =
           - { \xyz_inst.f07[1] , \xyz_inst.f06[1] , \xyz_inst.f05[1] , \xyz_inst.f04[1] ,
               \xyz_inst.f03[1] , 3'h0 };

endmodule
  )";

  std::istringstream in_ntk( file );

  bound_network ntk( gates );
  const auto result_ntk = rinox::io::verilog::read_verilog( in_ntk, rinox::io::reader( ntk ) );

  /* structural checks */
  CHECK( result_ntk == lorina::return_code::success );
  CHECK( ntk.num_pis() == 4 );
  CHECK( ntk.num_pos() == 8 );
  CHECK( ntk.size() == 20 ); // 2 constants, 4 PIs, 2 buffers
  CHECK( ntk.num_gates() == 14 );

  std::ostringstream out;
  rinox::io::verilog::write_verilog( ntk, out );

  std::string expected =
      "module top( a , b , xyz_inst.out1 );\n"
      "  input [0:1] a ;\n"
      "  input [0:1] b ;\n"
      "  output [0:7] xyz_inst.out1 ;\n"
      "  wire n0 , n1 , n6 , n7 , n8 , n9 , n10 , n11 ;\n"
      "  zero  g00( .O (n0) );\n"
      "  one   g01( .O (n1) );\n"
      "  inv2  g02( .a (n0), .O (xyz_inst.out1[0]) );\n"
      "  inv2  g03( .a (n0), .O (xyz_inst.out1[1]) );\n"
      "  inv2  g04( .a (n0), .O (xyz_inst.out1[2]) );\n"
      "  buf   g05( .a (b[0]), .O (n11) );\n"
      "  inv2  g06( .a (n11), .O (xyz_inst.out1[3]) );\n"
      "  buf   g07( .a (a[1]), .O (n10) );\n"
      "  inv2  g08( .a (n10), .O (xyz_inst.out1[4]) );\n"
      "  inv1  g09( .a (a[0]), .O (n8) );\n"
      "  xor2  g10( .a (n8), .b (b[1]), .O (n9) );\n"
      "  inv2  g11( .a (n9), .O (xyz_inst.out1[5]) );\n"
      "  xor2  g12( .a (a[1]), .b (b[1]), .O (n7) );\n"
      "  inv2  g13( .a (n7), .O (xyz_inst.out1[6]) );\n"
      "  and2  g14( .a (a[0]), .b (b[0]), .O (n6) );\n"
      "  inv2  g15( .a (n6), .O (xyz_inst.out1[7]) );\n"
      "endmodule\n";

  CHECK( out.str() == expected );
}

TEST_CASE( "Corner case 5 - signal names with erroneous spaces and vector instantiation of inputs and outputs", "[verilog_reader]" )
{

  using bound_network = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  std::vector<mockturtle::gate> gates;

  std::istringstream in_lib( test_library );
  auto result_lib = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result_lib == lorina::return_code::success );

  std::string file =
      R"(module top(
    \ripple_inst.a[0], \ripple_inst.a[1] ,
    \ripple_inst.o[0], \ripple_inst.o[1]
);
    // Port directions with vector widths
    input  \ripple_inst.a[0], \ripple_inst.a[1] ;
    output \ripple_inst.o[0], \ripple_inst.o[1] ;

    // Internal nets (all used)
    wire _008_;
    wire _038_;
    wire t0;

    // Your exact assignments (note the spaces before [1])
    assign _008_ = \ripple_inst.a [1];
    assign \ripple_inst.o [1] = _038_;
    
    // A couple of small gates to avoid floaters
    // Build _038_ from a bit of the input bus and _008_
    and2 g0 (.a(\ripple_inst.a [0]), .b(_008_), .O(_038_));
    // Drive the other output bit directly (also with spaces)
    buf  g1 (.a(\ripple_inst.a [0]), .O(\ripple_inst.o [0]));

endmodule)";

  std::istringstream in_ntk( file );

  bound_network ntk( gates );
  const auto result_ntk = rinox::io::verilog::read_verilog( in_ntk, rinox::io::reader( ntk ) );

  /* structural checks */
  CHECK( result_ntk == lorina::return_code::success );
  CHECK( ntk.num_pis() == 2 );
  CHECK( ntk.num_pos() == 2 );
  CHECK( ntk.size() == 6 ); // 2 constants, 4 PIs, 2 buffers
  CHECK( ntk.num_gates() == 2 );

  std::ostringstream out;
  rinox::io::verilog::write_verilog( ntk, out );

  std::string expected =
      "module top( ripple_inst.a , ripple_inst.o );\n"
      "  input [0:1] ripple_inst.a ;\n"
      "  output [0:1] ripple_inst.o ;\n"
      "  wire n0 , n1 ;\n"
      "  zero  g0( .O (n0) );\n"
      "  one   g1( .O (n1) );\n"
      "  buf   g2( .a (ripple_inst.a[0]), .O (ripple_inst.o[0]) );\n"
      "  and2  g3( .a (ripple_inst.a[0]), .b (ripple_inst.a[1]), .O (ripple_inst.o[1]) );\n"
      "endmodule\n";

  CHECK( out.str() == expected );
}