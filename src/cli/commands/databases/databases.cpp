#include "databases.hpp"
#include "../../context.hpp"
#include <fstream>
#include <iostream>
#include <lorina/genlib.hpp>
#include <rinox/databases/database_generator.hpp>
#include <rinox/databases/databases.hpp>
#include <rinox/libraries/augmented_library.hpp>

// Minimal build hook—replace with your actual builder.
static bool build_database( CLIContext& ctx,
                            rinox::databases::database_gen_params ps )
{
  if ( ps.num_vars == 4u )
  {
    static constexpr rinox::network::design_type_t design_t = rinox::network::design_type_t::CELL_BASED;
    rinox::databases::database_generator<design_t, 4u, 2u> db_gen( ctx.gates );
    db_gen.run( ps );
    std::cout << "Generated database" << std::endl;
    auto tmp = db_gen.extract_db();
    std::cout << "Generated database" << std::endl;
    ctx.db4.emplace( std::move( tmp ) );
    std::cout << "Assigned database" << std::endl;
    return true;
  }
  else
  {
    std::cerr << "Only num-vars = 4 for now.\n";
    return false;
  }
}

static void print_make_db_usage()
{
  std::cerr << "Usage: make_db --method <string> --num-vars <N> --metric <area|delay|power>\n"
               "Examples:\n"
               "  make_db --method mapp --num-vars 4 --metric area\n"
               "  make_db --method mapp --num-vars 4 --metric delay\n"
               "  make_db --method mapp --num-vars 4 --metric power\n"
               "Default:\n"
               "  make_db = (make_db --method mapp --num-vars 4 --metric area)\n";
}

static void cmd_make_database( CLIContext& ctx, const std::vector<std::string>& args )
{
  if ( args.empty() || std::find( args.begin(), args.end(), "--help" ) != args.end() )
  {
    print_make_db_usage();
    return;
  }

  if ( ctx.gates.empty() )
  {
    std::cerr << "Error: load a library first with `read_genlib <file.genlib>`.\n";
    return;
  }

  rinox::databases::database_gen_params dbps;

  // simple argv-style parsing
  for ( size_t i = 0; i < args.size(); ++i )
  {
    const std::string& a = args[i];
    auto need_val = [&]( const char* flag ) -> bool {
      if ( i + 1 >= args.size() )
      {
        std::cerr << "Error: " << flag << " requires a value.\n";
        return false;
      }
      return true;
    };

    if ( a == "--method" )
    {
      if ( !need_val( "--method" ) )
        return;
      dbps.method = args[++i];
    }
    else if ( a == "--num-vars" )
    {
      if ( !need_val( "--num-vars" ) )
        return;
      try
      {
        long long v = std::stoll( args[++i] );
        if ( v < 1 || v > 64 )
        { // adjust bounds as you like
          std::cerr << "Error: --num-vars must be in [1, 64].\n";
          return;
        }
        dbps.num_vars = static_cast<uint32_t>( v );
      }
      catch ( ... )
      {
        std::cerr << "Error: --num-vars expects an integer.\n";
        return;
      }
    }
    else if ( a == "--metric" )
    {
      if ( !need_val( "--metric" ) )
        return;
      dbps.metric = args[++i];
    }
    else if ( a.rfind( "--", 0 ) == 0 )
    {
      std::cerr << "Warning: unknown option '" << a << "' (ignored).\n";
    }
    else
    {
      std::cerr << "Warning: stray argument '" << a << "' (ignored).\n";
    }
  }

  // build
  if ( !build_database( ctx, dbps ) )
  {
    std::cerr << "Error: database build failed.\n";
    return;
  }

  std::cout << "Database built: method=" << dbps.method
            << " num_vars=" << dbps.num_vars
            << " metric=" << dbps.metric << "\n";
}

static void cmd_dump_database( CLIContext& ctx, const std::vector<std::string>& args )
{
  if ( args.size() < 2 )
  {
    std::cerr << "Usage: dump_db <filename>.v\n";
    return;
  }
  if ( !ctx.db4 )
  {
    std::cerr << "No database generated\n";
    return;
  }
  ( *( ctx.db4 ) ).commit( args[1] );
}

static void cmd_read_database( CLIContext& ctx, const std::vector<std::string>& args )
{
  if ( args.size() < 2 )
  {
    std::cerr << "Usage: read_db <filename>.v\n";
    return;
  }

  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using Db = rinox::databases::mapped_database<Ntk, 4u>;
  rinox::libraries::augmented_library<rinox::network::design_type_t::CELL_BASED> lib( ctx.gates );
  Db tmp( lib );
  tmp.load( args[1] );
  ctx.db4.emplace( std::move( tmp ) );
}

std::map<std::string, CommandHandler> register_db_commands()
{
  return {
      { "make_db", cmd_make_database },
      { "dump_db", cmd_dump_database },
      { "read_db", cmd_read_database } };
}
