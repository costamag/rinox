#pragma once
#include <cstdio>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <stdexcept>
#include <filesystem>

namespace rinox::experiments::array
{

template<class Params>
struct area_resyn_t
{
  Params ps;
  std::string input;

  struct {
    std::vector<std::string> todo;
    std::filesystem::path path;
  } database;
};

template<class Params>
inline area_resyn_t<Params> load_area_resyn_params( rapidjson::Value const& doc, lorina::diagnostic_engine* diag )
{
  area_resyn_t<Params> res;
  Params& ps = res.ps;

  // Ensure "resynthesis" exists
  if ( !doc.HasMember( "resyn" ) || !doc["resyn"].IsObject() )
    return res;

  const auto& rp = doc["resyn"];

  // -----------------------------
  // Parse "profiler_ps"
  // -----------------------------
  if ( rp.HasMember( "kStopAfterNth" ) && rp["kStopAfterNth"].IsUint() )
  {
    ps.kStopAfterNth = rp["kStopAfterNth"].GetUint();
  }

  return res;
}

} // namespace rinox::experiments
