#pragma once
#include <cstdio>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <stdexcept>
#include <filesystem>

namespace rinox::experiments
{

template<uint32_t MaxNumLeaves>
struct windowing_t
{
  static constexpr uint32_t max_num_leaves = MaxNumLeaves;
  bool preserve_depth = false;
  int odc_levels = 0;
  uint32_t skip_fanout_limit_for_divisors = 200;
  uint32_t max_num_divisors = 1000;
  uint32_t max_reconvergence = 10;
};

template<uint32_t MaxNumLeaves>
inline windowing_t<MaxNumLeaves> load_windowing_params( rapidjson::Value const& doc, lorina::diagnostic_engine* diag )
{
  windowing_t<MaxNumLeaves> res;
  if ( !doc.HasMember( "windowing" ) || !doc["windowing"].IsObject() )
    return res;

  const auto& wm = doc["windowing"];

  // -----------------------------
  // Parse "window_manager_ps"
  // -----------------------------
  if ( wm.HasMember( "preserve_depth" ) && wm["preserve_depth"].IsBool() )
    res.preserve_depth = wm["preserve_depth"].GetBool();
  if ( wm.HasMember( "odc_levels" ) && wm["odc_levels"].IsInt() )
    res.odc_levels = wm["odc_levels"].GetInt();
  if ( wm.HasMember( "skip_fanout_limit_for_divisors" ) && wm["skip_fanout_limit_for_divisors"].IsUint() )
    res.skip_fanout_limit_for_divisors = wm["skip_fanout_limit_for_divisors"].GetUint();
  if ( wm.HasMember( "max_num_divisors" ) && wm["max_num_divisors"].IsUint() )
    res.max_num_divisors = wm["max_num_divisors"].GetUint();
  if ( wm.HasMember( "max_reconvergence" ) && wm["max_reconvergence"].IsUint() )
    res.max_reconvergence = wm["max_reconvergence"].GetUint();
  return res;
}

} // namespace rinox::experiments
