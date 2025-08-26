#include "tech_map.hpp"
#include <cstdio>
#include <rinox/diagnostics.hpp>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <mockturtle/algorithms/emap.hpp>
#include <stdexcept>

namespace rinox::experiments
{

tech_map_t load_tech_map_params( rapidjson::Document const& doc, lorina::diagnostic_engine * diag )
{
  tech_map_t cfg;

  if (!doc.HasMember("tech-map") || !doc["tech-map"].IsObject())
  {
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::fatal, "`tech-map` section missing or invalid");
    return cfg;
  }

  const auto& tech_map = doc["tech-map"];

  // -- Parse library params --
  if (tech_map.HasMember("lib-params") && tech_map["lib-params"].IsObject())
  {
    const auto& lib = tech_map["lib-params"];

    if (lib.HasMember("ignore_symmetries") && lib["ignore_symmetries"].IsBool())
      cfg.tps.ignore_symmetries = lib["ignore_symmetries"].GetBool();

    if (lib.HasMember("verbose") && lib["verbose"].IsBool())
      cfg.tps.verbose = lib["verbose"].GetBool();
  }

  // -- Parse mapping params --
  if (tech_map.HasMember("map-params") && tech_map["map-params"].IsObject())
  {
    const auto& map = tech_map["map-params"];

    if (map.HasMember("cut_limit") && map["cut_limit"].IsUint())
      cfg.mps.cut_enumeration_ps.cut_limit = map["cut_limit"].GetUint();

    if (map.HasMember("minimize_truth_table") && map["minimize_truth_table"].IsBool())
      cfg.mps.cut_enumeration_ps.minimize_truth_table = map["minimize_truth_table"].GetBool();

    if (map.HasMember("matching_mode") && map["matching_mode"].IsString())
    {
      std::string mode = map["matching_mode"].GetString();

      if (mode == "boolean")
        cfg.mps.matching_mode = mockturtle::emap_params::matching_mode_t::boolean;
      else if (mode == "structural")
        cfg.mps.matching_mode = mockturtle::emap_params::matching_mode_t::structural;
      else if (mode == "hybrid")
        cfg.mps.matching_mode = mockturtle::emap_params::matching_mode_t::hybrid;
      else
        rinox::diagnostics::REPORT_DIAG(diag, lorina::diagnostic_level::warning,
                                        "Unknown `matching_mode`: {}", mode);
    }

    if (map.HasMember("area_oriented_mapping") && map["area_oriented_mapping"].IsBool())
      cfg.mps.area_oriented_mapping = map["area_oriented_mapping"].GetBool();

    if (map.HasMember("map_multioutput") && map["map_multioutput"].IsBool())
      cfg.mps.map_multioutput = map["map_multioutput"].GetBool();

    if (map.HasMember("relax_required") && map["relax_required"].IsUint())
      cfg.mps.relax_required = map["relax_required"].GetUint();

    if (map.HasMember("required_time") && map["required_time"].IsUint())
      cfg.mps.required_time = map["required_time"].GetUint();

    if (map.HasMember("area_flow_rounds") && map["area_flow_rounds"].IsUint())
      cfg.mps.area_flow_rounds = map["area_flow_rounds"].GetUint();

    if (map.HasMember("ela_rounds") && map["ela_rounds"].IsUint())
      cfg.mps.ela_rounds = map["ela_rounds"].GetUint();

    if (map.HasMember("eswp_rounds") && map["eswp_rounds"].IsUint())
      cfg.mps.eswp_rounds = map["eswp_rounds"].GetUint();

    if (map.HasMember("switching_activity_patterns") && map["switching_activity_patterns"].IsUint())
      cfg.mps.switching_activity_patterns = map["switching_activity_patterns"].GetUint();

    if (map.HasMember("use_match_alternatives") && map["use_match_alternatives"].IsBool())
      cfg.mps.use_match_alternatives = map["use_match_alternatives"].GetBool();

    if (map.HasMember("remove_dominated_cuts") && map["remove_dominated_cuts"].IsBool())
      cfg.mps.remove_dominated_cuts = map["remove_dominated_cuts"].GetBool();

    if (map.HasMember("remove_overlapping_multicuts") && map["remove_overlapping_multicuts"].IsBool())
      cfg.mps.remove_overlapping_multicuts = map["remove_overlapping_multicuts"].GetBool();

    if (map.HasMember("verbose") && map["verbose"].IsBool())
      cfg.mps.verbose = map["verbose"].GetBool();
  }
  else
  {
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::fatal, "`map-params` missing from `tech-map`");
  }

  return cfg;
}

} // namespace rinox::experiments
