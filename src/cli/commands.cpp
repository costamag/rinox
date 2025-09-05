#include "commands.hpp"
#include "commands/databases/databases.hpp"
#include "commands/io/io.hpp"
#include "commands/libraries.hpp"
#include "commands/opto/resynthesis.hpp"

std::map<std::string, CommandHandler> register_commands()
{
  auto cmds = register_library_commands();
  auto io_cmds = register_io_commands();
  cmds.insert( io_cmds.begin(), io_cmds.end() );
  auto db_cmds = register_db_commands();
  cmds.insert( db_cmds.begin(), db_cmds.end() );
  /* opto */
  auto opto_cmds = register_resyn_commands();
  cmds.insert( opto_cmds.begin(), opto_cmds.end() );
  return cmds;
}
