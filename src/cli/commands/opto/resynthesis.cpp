#include "resynthesis.hpp"
#include "../../../../include/rinox/opto/algorithms/resynthesize.hpp"
#include "../../../../include/rinox/network/network.hpp"
#include "../../context.hpp"
#include <fstream>
#include <iostream>
#include <lorina/verilog.hpp>

static void cmd_resyn( CLIContext& ctx, const std::vector<std::string>& args )
{
  if ( !ctx.ntk )
  {
    std::cerr << "Error: load a network first with read_verilog.\n";
    return;
  }
  std::string target = "area";

  
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::CELL_BASED, 2>;
  using Db = rinox::databases::mapped_database<Ntk, 4u>;
  using Ps = rinox::opto::algorithms::default_resynthesis_params<12>;
  Ps ps;

  for ( size_t i = 1; i < args.size(); ++i )
  {
    const auto& arg = args[i];

    if ( arg == "--target" && i + 1 < args.size() )
    {
      target = std::string(args[++i]);
    }
    else if ( arg == "--odc-levels" && i + 1 < args.size() )
    {
      ps.window_manager_ps.odc_levels = std::stoi( args[++i] );
    }
    else if ( arg == "--num-divs" && i + 1 < args.size() )
    {
      ps.window_manager_ps.max_num_divisors = std::stoi( args[++i] );
    }
    else if ( arg == "--stop-at" && i + 1 < args.size() )
    {
      ps.stop_after_nth = std::stoul( args[++i] );
    }
    else if ( arg == "--rewire" )
    {
      ps.try_rewire = true;
    }
    else if ( arg == "--window" )
    {
      ps.try_window = true;
    }
    else if ( arg == "--struct" )
    {
      ps.try_struct = true;
    }
    else if ( arg == "--simula" )
    {
      ps.try_simula = true;
    }
    else if ( arg == "--help" )
    {
      std::cout << "Usage: resyn [options]\n\n"
                << "Options:\n"
                << "  --target {area,delay,power}   Optimization target (default: area)\n"
                << "  --odc-levels N                ODC levels [0–2] (default: 0)\n"
                << "  --num-divs N                  Max number of divisors (default: 100)\n"
                << "  --stop-at N                   Stop after N successful rewrites\n"
                << "  --rewire                      Enable fanin rewiring\n"
                << "  --window                      Enable window-based rewriting\n"
                << "  --struct                      Enable structural cut rewriting\n"
                << "  --simula                      Enable simulation-guided rewriting\n";
      return;
    }
    else
    {
      std::cerr << "Unknown option: " << arg << "\n";
    }
  }


  if ( !ctx.db4 )
  {
    std::cerr << "No database loaded\n";
    return;
  }


  rinox::opto::algorithms::resynthesis_stats st;
  if ( target == "area" )
  {
    if ( ctx.db4 )
    {
      rinox::opto::algorithms::resynthesis_stats st;
      rinox::opto::algorithms::area_resynthesize<Ntk, Db, Ps>( *ctx.ntk, *ctx.db4, ps, &st );
    }
    return;
  }
  else if ( target == "delay" )
  {
    if ( ctx.db4 )
    {
      rinox::opto::algorithms::resynthesis_stats st;
      rinox::opto::algorithms::delay_resynthesize<Ntk, Db, Ps>( *ctx.ntk, *ctx.db4, ps, &st );
    }
    return;
  }
  else if ( target == "power" )
  {
    if ( ctx.db4 )
    {
      rinox::opto::algorithms::resynthesis_stats st;
      rinox::opto::algorithms::power_resynthesize<Ntk, Db, Ps>( *ctx.ntk, *ctx.db4, ps, &st );
    }
    return;
  }
  else
  {
    std::cerr << "Invalid target: " << target << "\n";
  }
}


std::map<std::string, CommandHandler> register_resyn_commands()
{
  return {
      { "resyn", cmd_resyn }
  };
}
