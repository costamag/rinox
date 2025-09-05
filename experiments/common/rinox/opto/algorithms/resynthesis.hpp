#pragma once
#include <cstdio>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <stdexcept>
#include <filesystem>

namespace rinox::experiments
{

template<class Params>
struct resynthesis_t
{
  Params ps;
  std::string input;

  struct {
    std::vector<std::string> todo;
    std::string tech;
    std::filesystem::path path;
  } database;
};

template<class Params>
inline resynthesis_t<Params> load_resynthesis_params( rapidjson::Value const& doc, lorina::diagnostic_engine* diag )
{
  resynthesis_t<Params> res;
  Params& ps = res.ps;

  // Ensure "resynthesis" exists
  if ( !doc.HasMember( "resynthesis" ) || !doc["resynthesis"].IsObject() )
    return res;

  const auto& rp = doc["resynthesis"];

  // -----------------------------
  // Parse "profiler_ps"
  // -----------------------------
  if ( rp.HasMember( "profiler_ps" ) && rp["profiler_ps"].IsObject() )
  {
    const auto& p = rp["profiler_ps"];
    if ( p.HasMember( "max_num_roots" ) && p["max_num_roots"].IsUint() )
      ps.profiler_ps.max_num_roots = p["max_num_roots"].GetUint();

    if ( p.HasMember( "input_arrivals" ) && p["input_arrivals"].IsDouble() )
      ps.profiler_ps.input_arrivals = { p["input_arrivals"].GetDouble() };

    if ( p.HasMember( "output_required" ) )
    {
      if ( p["output_required"].IsString() && std::string( p["output_required"].GetString() ) == "INF" )
        ps.profiler_ps.output_required = { -1 };
      else if ( p["output_required"].IsDouble() )
        ps.profiler_ps.output_required = { p["output_required"].GetDouble() };
    }

    if ( p.HasMember( "eps" ) && p["eps"].IsNumber() )
      ps.profiler_ps.eps = p["eps"].GetDouble();
  }

  // -----------------------------
  // Parse "window_manager_ps"
  // -----------------------------
  if ( rp.HasMember( "window_manager_ps" ) && rp["window_manager_ps"].IsObject() )
  {
    const auto& wm = rp["window_manager_ps"];
    if ( wm.HasMember( "preserve_depth" ) && wm["preserve_depth"].IsBool() )
      ps.window_manager_ps.preserve_depth = wm["preserve_depth"].GetBool();
    if ( wm.HasMember( "odc_levels" ) && wm["odc_levels"].IsInt() )
      ps.window_manager_ps.odc_levels = wm["odc_levels"].GetInt();
    if ( wm.HasMember( "skip_fanout_limit_for_divisors" ) && wm["skip_fanout_limit_for_divisors"].IsUint() )
      ps.window_manager_ps.skip_fanout_limit_for_divisors = wm["skip_fanout_limit_for_divisors"].GetUint();
    if ( wm.HasMember( "max_num_divisors" ) && wm["max_num_divisors"].IsUint() )
      ps.window_manager_ps.max_num_divisors = wm["max_num_divisors"].GetUint();
  }

  // -----------------------------
  // Parse other flags
  // -----------------------------
  if ( rp.HasMember("params") && rp["params"].IsObject() )
  {
    const auto& params = rp["params"];

    if ( params.HasMember( "use_dont_cares" ) && params["use_dont_cares"].IsBool() )
      ps.use_dont_cares = params["use_dont_cares"].GetBool();
    if ( params.HasMember( "try_rewire" ) && params["try_rewire"].IsBool() )
      ps.try_rewire = params["try_rewire"].GetBool();
    if ( params.HasMember( "try_struct" ) && params["try_struct"].IsBool() )
      ps.try_struct = params["try_struct"].GetBool();
    if ( params.HasMember( "try_window" ) && params["try_window"].IsBool() )
      ps.try_window = params["try_window"].GetBool();
    if ( params.HasMember( "try_simula" ) && params["try_simula"].IsBool() )
      ps.try_simula = params["try_simula"].GetBool();
    if ( params.HasMember( "dynamic_database" ) && params["dynamic_database"].IsBool() )
      ps.dynamic_database = params["dynamic_database"].GetBool();
    if ( params.HasMember( "fanout_limit" ) && params["fanout_limit"].IsUint() )
      ps.fanout_limit = params["fanout_limit"].GetUint();
    if ( params.HasMember( "stop_after_nth" ) && params["stop_after_nth"].IsUint() )
      ps.stop_after_nth = params["stop_after_nth"].GetUint();
    // -----------------------------
    // Parse "database" section
    // -----------------------------
    if ( rp.HasMember( "database" ) && rp["database"].IsObject() )
    {
      const auto& db = rp["database"];
  
      // Parse "todo"
      if ( db.HasMember("todo") && db["todo"].IsArray() )
      {
        for ( const auto& entry : db["todo"].GetArray() )
        {
          if ( entry.IsString() )
            res.database.todo.emplace_back( entry.GetString() );
          else
            rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::warning,
              "Invalid entry in `database.todo` (must be a string)");
        }
      }
  
      if ( res.database.todo[0] == "load" )
      {
        // Parse "type"
        if ( db.HasMember("tech") && db["tech"].IsString() )
          res.database.tech = db["tech"].GetString();
        else
          rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::warning,
            "`database.tech` missing or invalid");
  
        std::string name;
        // Parse "name"
        if ( db.HasMember("name") && db["name"].IsString() )
          name = db["name"].GetString();
        else
          rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::warning,
            "`database.name` missing or invalid");
  
        std::filesystem::path databases = fs::path(RINOX_SOURCE_DIR) / "databases";
        std::string filename = name + ".v";
        res.database.path = databases / res.database.tech / filename;  
      }
      else
      {
        rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::fatal,
          "`database` section missing from `resynthesis`");
      }
    }
  }

  return res;
}

} // namespace rinox::experiments
