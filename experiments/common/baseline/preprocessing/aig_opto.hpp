#pragma once
#include <rinox/traits.hpp>
#include <rinox/diagnostics/diagnostics.hpp>
#include <mockturtle/networks/aig.hpp>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <filesystem>
#include <set>
#include <string>
#include <vector>

namespace rinox::experiments
{

struct aig_opto_t {
  std::string name;
  std::string algorithm;
  int iterations = 1;
  rapidjson::Value const* params; // raw pointer into parsed JSON
};

std::vector<aig_opto_t> parse_aig_opto( rapidjson::Document const& doc, lorina::diagnostic_engine* diag = nullptr );

void apply_aig_opto(mockturtle::aig_network&, const aig_opto_t& );

} // namespace rinox::experiments
