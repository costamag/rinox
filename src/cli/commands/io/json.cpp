#include "json.hpp"
#include "../../../../include/rinox/io/json/json.hpp"
#include <rinox/diagnostics/diagnostics.hpp>
#include "../../../../include/rinox/network/network.hpp"
#include "../../context.hpp"
#include <fstream>
#include <iostream>

static void cmd_read_json( CLIContext& ctx, const std::vector<std::string>& args )
{
  if ( args.size() < 2 )
  {
    std::cerr << "Usage: read_json <filename>\n";
    return;
  }
  if ( ctx.gates.empty() )
  {
    std::cerr << "Error: load a library first with read_genlib.\n";
    return;
  }
  std::ifstream in_ntk( args[1] );
  if ( !in_ntk )
  {
    std::cerr << "Cannot open " << args[1] << "\n";
    return;
  }
  ctx.ntk.emplace( ctx.gates );

  rinox::diagnostics::text_diagnostics consumer;
  lorina::diagnostic_engine diag( &consumer );

  auto result_ntk = rinox::io::json::read_json(
      in_ntk, rinox::io::reader( *ctx.ntk ), &diag );
  if ( result_ntk == lorina::return_code::success )
    std::cout << "Design loaded.\n";
  else
    std::cerr << "Failed to read json.\n";
}

static void cmd_write_json( CLIContext& ctx, const std::vector<std::string>& args )
{
  if ( !ctx.ntk )
  {
    std::cerr << "Error: load a network first with read_verilog or read_json.\n";
    return;
  }
  if ( args.size() < 2 )
  {
    std::cerr << "Usage: write_json <filename>\n";
    return;
  }
  std::ofstream out( args[1] );
  if ( !out )
  {
    std::cerr << "Cannot write to " << args[1] << "\n";
    return;
  }
  rinox::io::json::write_json( *ctx.ntk, out );
  std::cout << "Design written to " << args[1] << "\n";
}

std::map<std::string, CommandHandler> register_json_commands()
{
  return {
      { "read_json", cmd_read_json },
      { "write_json", cmd_write_json } };
}
