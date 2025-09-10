#pragma once
#include <cstdio>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <stdexcept>
#include <lorina/lorina.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>

namespace rinox::experiments
{

struct lut_map_t
{
  mockturtle::lut_map_params mps;
};

lut_map_t load_lut_map_params( rapidjson::Document const& doc, lorina::diagnostic_engine * diag );

} // namespace rinox::experiments
