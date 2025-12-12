#pragma once
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include <cstdio>
#include <filesystem>
#include <stdexcept>

namespace rinox::experiments::array {

template <class Params>
struct area_struct_t {
  Params ps;
  std::string input;

  struct {
    std::vector<std::string> todo;
    std::filesystem::path path;
  } database;
};

template <class Params>
inline area_struct_t<Params> load_area_struct_params(
    rapidjson::Value const& doc, lorina::diagnostic_engine* diag) {
  area_struct_t<Params> res;
  Params& ps = res.ps;

  // Ensure "resynthesis" exists
  if (!doc.HasMember("resyn") || !doc["resyn"].IsObject()) return res;

  const auto& rp = doc["resyn"];

  // -----------------------------
  // Parse "profiler_ps"
  // -----------------------------
  if (rp.HasMember("stop_after_nth") && rp["stop_after_nth"].IsUint()) {
    ps.stop_after_nth = rp["stop_after_nth"].GetUint();
  }
  if (rp.HasMember("max_num_divisors") && rp["max_num_divisors"].IsUint()) {
    ps.max_num_divisors = rp["max_num_divisors"].GetUint();
  }

  if (rp.HasMember("max_reconvergence") && rp["max_reconvergence"].IsUint()) {
    ps.max_reconvergence = rp["max_reconvergence"].GetUint();
  }
  if (rp.HasMember("num_levels_tfo") && rp["num_levels_tfo"].IsUint()) {
    ps.max_num_levels_tfo = rp["num_levels_tfo"].GetUint();
  }

  return res;
}

}  // namespace rinox::experiments::array
