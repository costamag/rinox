#include <mockturtle/io/genlib_reader.hpp>
#include <rinox/databases/mapped_database.hpp>
#include <rinox/databases/database_generator.hpp>
#include "../../../common/rinox/databases/mapped_database.hpp"

int main()
{
  /* parse the configuration file */
  rinox::diagnostics::text_diagnostics client;
  lorina::diagnostic_engine diag( &client );
  std::string const exp_path = std::string( RINOX_EXPERIMENTS_DIR ) + "/studies/databases/map_pclasses/";
  std::string const config_path = exp_path + "config.json";
  auto doc = rinox::experiments::load_json_doc( config_path, &diag );
  if ( !doc ) return 1;
  
  /* parse the configuration file */
  auto const& dbc = rinox::experiments::load_database_gen( *doc, &diag );

  /* load the gtechnology library */
  fmt::print( "[i] processing technology library\n" );
  std::vector<mockturtle::gate> gates;
  std::ifstream in( dbc.technology.path );
  if ( lorina::read_genlib( in, mockturtle::genlib_reader( gates ) ) != lorina::return_code::success )
    return 1;

  /* generate the database */
  rinox::databases::database_generator<rinox::network::design_type_t::CELL_BASED, 4u, 2u> gen( gates );
  /* load the database */
  rinox::databases::database_gen_params ps;
  ps.method = dbc.parameters.method;
  ps.num_vars = dbc.parameters.num_vars;
  ps.multiple_candidates = dbc.parameters.multiple_candidates;
  ps.verbose = dbc.parameters.verbose;
  ps.metric = dbc.parameters.metric;
  ps.output_path = dbc.path_out + dbc.technology.name + "_" + ps.method + "_" + std::to_string( ps.num_vars ) + ".v";
  std::cout << ps.output_path << std::endl;
  gen.run( ps );

  return 0;
}