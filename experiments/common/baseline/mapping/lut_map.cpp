#include "lut_map.hpp"
#include <cstdio>
#include <rinox/diagnostics/diagnostics.hpp>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <stdexcept>

namespace rinox::experiments
{

lut_map_t load_lut_map_params( rapidjson::Document const& doc, lorina::diagnostic_engine * diag )
{
  lut_map_t cfg;

  if ( !doc.HasMember("lut-map") || !doc["lut-map"].IsObject() )
  {
    rinox::diagnostics::REPORT_DIAG_RAW(
      diag, lorina::diagnostic_level::fatal,
      "`lut-map` section missing or invalid" );
    return cfg;
  }

  const auto& lut_map = doc["lut-map"];

  // -- Parse mapping params --
  if ( lut_map.HasMember("map-params") && lut_map["map-params"].IsObject() )
  {
    const auto& map = lut_map["map-params"];

    if ( map.HasMember("cut_size") && map["cut_size"].IsUint() )
      cfg.mps.cut_enumeration_ps.cut_size = map["cut_size"].GetUint();

    if ( map.HasMember("cut_limit") && map["cut_limit"].IsUint() )
      cfg.mps.cut_enumeration_ps.cut_limit = map["cut_limit"].GetUint();

    if ( map.HasMember("minimize_truth_table") && map["minimize_truth_table"].IsBool() )
      cfg.mps.cut_enumeration_ps.minimize_truth_table = map["minimize_truth_table"].GetBool();

    if ( map.HasMember("area_oriented_mapping") && map["area_oriented_mapping"].IsBool() )
      cfg.mps.area_oriented_mapping = map["area_oriented_mapping"].GetBool();

    if ( map.HasMember("required_delay") && map["required_delay"].IsUint() )
      cfg.mps.required_delay = map["required_delay"].GetUint();

    if ( map.HasMember("relax_required") && map["relax_required"].IsUint() )
      cfg.mps.relax_required = map["relax_required"].GetUint();

    if ( map.HasMember("recompute_cuts") && map["recompute_cuts"].IsBool() )
      cfg.mps.recompute_cuts = map["recompute_cuts"].GetBool();

    if ( map.HasMember("area_share_rounds") && map["area_share_rounds"].IsUint() )
      cfg.mps.area_share_rounds = map["area_share_rounds"].GetUint();

    if ( map.HasMember("area_flow_rounds") && map["area_flow_rounds"].IsUint() )
      cfg.mps.area_flow_rounds = map["area_flow_rounds"].GetUint();

    if ( map.HasMember("ela_rounds") && map["ela_rounds"].IsUint() )
      cfg.mps.ela_rounds = map["ela_rounds"].GetUint();

    if ( map.HasMember("edge_optimization") && map["edge_optimization"].IsBool() )
      cfg.mps.edge_optimization = map["edge_optimization"].GetBool();

    if ( map.HasMember("cut_expansion") && map["cut_expansion"].IsBool() )
      cfg.mps.cut_expansion = map["cut_expansion"].GetBool();

    if ( map.HasMember("remove_dominated_cuts") && map["remove_dominated_cuts"].IsBool() )
      cfg.mps.remove_dominated_cuts = map["remove_dominated_cuts"].GetBool();

    if ( map.HasMember("collapse_mffcs") && map["collapse_mffcs"].IsBool() )
      cfg.mps.collapse_mffcs = map["collapse_mffcs"].GetBool();

    if ( map.HasMember("sop_balancing") && map["sop_balancing"].IsBool() )
      cfg.mps.sop_balancing = map["sop_balancing"].GetBool();

    if ( map.HasMember("esop_balancing") && map["esop_balancing"].IsBool() )
      cfg.mps.esop_balancing = map["esop_balancing"].GetBool();

    if ( map.HasMember("cost_cache_vars") && map["cost_cache_vars"].IsUint() )
      cfg.mps.cost_cache_vars = map["cost_cache_vars"].GetUint();

    if ( map.HasMember("verbose") && map["verbose"].IsBool() )
      cfg.mps.verbose = map["verbose"].GetBool();
  }
  else
  {
    rinox::diagnostics::REPORT_DIAG_RAW(
      diag, lorina::diagnostic_level::fatal,
      "`map-params` missing from `lut-map`" );
  }

  return cfg;
}

} // namespace rinox::experiments
