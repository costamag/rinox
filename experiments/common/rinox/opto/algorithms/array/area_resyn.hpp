#pragma once
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include <cstdio>
#include <filesystem>
#include <stdexcept>

namespace rinox::experiments::array {

template <class Params>
struct area_resyn_t {
  Params ps;
  std::string input;

  struct {
    std::vector<std::string> todo;
    std::filesystem::path path;
  } database;
};

template <class Params>
inline area_resyn_t<Params> load_area_resyn_params(
    rapidjson::Value const& doc, lorina::diagnostic_engine* diag) {
  area_resyn_t<Params> res;
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
  if (rp.HasMember("odc_levels") && rp["odc_levels"].IsUint()) {
    ps.odc_levels = rp["odc_levels"].GetUint();
  }
  if (rp.HasMember("cec") && rp["cec"].IsBool()) {
    ps.cec = rp["cec"].GetBool();
  }
  if (rp.HasMember("try_0resub_solo") && rp["try_0resub_solo"].IsBool()) {
    ps.try_0resub_solo = rp["try_0resub_solo"].GetBool();
  }
  if (rp.HasMember("try_1resub_solo") && rp["try_1resub_solo"].IsBool()) {
    ps.try_1resub_solo = rp["try_1resub_solo"].GetBool();
  }
  if (rp.HasMember("try_2resub_solo") && rp["try_2resub_solo"].IsBool()) {
    ps.try_2resub_solo = rp["try_2resub_solo"].GetBool();
  }
  if (rp.HasMember("try_1resub_dual") && rp["try_1resub_dual"].IsBool()) {
    ps.try_1resub_dual = rp["try_1resub_dual"].GetBool();
  }
  if (rp.HasMember("try_1rexub_solo") && rp["try_1rexub_solo"].IsBool()) {
    ps.try_1rexub_solo = rp["try_1rexub_solo"].GetBool();
  }
  if (rp.HasMember("try_0rexub_solo") && rp["try_0rexub_solo"].IsBool()) {
    ps.try_0rexub_solo = rp["try_0rexub_solo"].GetBool();
  }
  if (rp.HasMember("try_merger_solo") && rp["try_merger_solo"].IsBool()) {
    ps.try_merger_solo = rp["try_merger_solo"].GetBool();
  }

  return res;
}

}  // namespace rinox::experiments::array
