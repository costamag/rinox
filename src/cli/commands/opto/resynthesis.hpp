#pragma once
#include "../../commands.hpp"
#include <map>
#include <string>

void resyn();

std::map<std::string, CommandHandler> register_resyn_commands();
