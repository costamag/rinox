#pragma once
#include "../../env.hpp"
#include <rinox/traits.hpp>
#include <rinox/diagnostics.hpp>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <mockturtle/io/genlib_reader.hpp>
#include <filesystem>
#include <set>
#include <string>
#include <vector>

namespace rinox::experiments
{

  struct database_gen_t
  {
    struct {
      std::string extension;
      std::string name;
      std::filesystem::path path;
  } technology;
  
  struct {
      std::string method;
      unsigned num_vars = 4;
      bool multiple_candidates = false;
      bool verbose = false;
      std::string metric;
  } parameters;

  std::filesystem::path path_out;
};

inline database_gen_t load_database_gen( rapidjson::Document const& doc,
                                          lorina::diagnostic_engine* diag )
{
  database_gen_t cfg;

  // ----------------------
  // Parse "technology"
  // ----------------------
  if ( doc.HasMember("technology") && doc["technology"].IsObject() )
  {
    const auto& t = doc["technology"];

    if ( t.HasMember("extension") && t["extension"].IsString() )
      cfg.technology.extension = t["extension"].GetString();
    else
      rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::warning,
        "`technology.extension` missing or invalid");

    if ( t.HasMember("name") && t["name"].IsString() )
      cfg.technology.name = t["name"].GetString();
    else
      rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::warning,
        "`technology.name` missing or invalid");
    
    cfg.path_out = std::filesystem::path( RINOX_SOURCE_DIR ) / "databases" / cfg.technology.name;
    std::string const file_name = cfg.technology.name + "." + cfg.technology.extension;
    cfg.technology.path = std::filesystem::path( RINOX_SOURCE_DIR ) / "techlib" / cfg.technology.extension / file_name;
  }
  else
  {
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::fatal,
      "`technology` section missing or invalid");
  }

  // ----------------------
  // Parse "parameters"
  // ----------------------
  if ( doc.HasMember("parameters") && doc["parameters"].IsObject() )
  {
    const auto& p = doc["parameters"];

    if ( p.HasMember("method") && p["method"].IsString() )
      cfg.parameters.method = p["method"].GetString();

    if ( p.HasMember("num_vars") && p["num_vars"].IsUint() )
      cfg.parameters.num_vars = p["num_vars"].GetUint();

    if ( p.HasMember("multiple_candidates") && p["multiple_candidates"].IsBool() )
      cfg.parameters.multiple_candidates = p["multiple_candidates"].GetBool();

    if ( p.HasMember("verbose") && p["verbose"].IsBool() )
      cfg.parameters.verbose = p["verbose"].GetBool();

    if ( p.HasMember("metric") && p["metric"].IsString() )
      cfg.parameters.metric = p["metric"].GetString();
  }
  else
  {
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::fatal,
      "`parameters` section missing or invalid");
  }

  return cfg;
}


} // namespace rinox::experiments
