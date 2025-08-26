#include "aig_opto.hpp"
#include <rinox/traits.hpp>
#include <rinox/diagnostics/diagnostics.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <filesystem>
#include <set>
#include <string>
#include <vector>

namespace rinox::experiments
{

std::vector<aig_opto_t> parse_aig_opto( rapidjson::Document const& doc, lorina::diagnostic_engine* diag )
{
  std::vector<aig_opto_t> pipeline;

  if (doc.HasMember("aig-opto") && doc["aig-opto"].IsObject())
  {  
    const auto& jb = doc["aig-opto"];
    for (auto it = jb.MemberBegin(); it != jb.MemberEnd(); ++it)
    {
      const std::string step_name = it->name.GetString();
      const auto& step_conf = it->value;

      if (!step_conf.HasMember("algorithm"))
        throw std::runtime_error("Missing algorithm field in step " + step_name);

      aig_opto_t step;
      step.name = step_name;
      step.algorithm = step_conf["algorithm"].GetString();
      step.iterations = step_conf.HasMember("iterations") ? step_conf["iterations"].GetInt() : 1;

      if (step_conf.HasMember("params") && step_conf["params"].IsObject())
        step.params = &step_conf["params"];
      else
        step.params = nullptr;

      pipeline.push_back(step);
    }
  }
  else {
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::warning, "no aig-opto");
    return pipeline;
  }

  return pipeline;
}

void apply_aig_balancing(mockturtle::aig_network& aig, rapidjson::Value const& jb)
{
  mockturtle::aig_balancing_params ps;

  if (jb.HasMember("minimize_levels") && jb["minimize_levels"].IsBool())
    ps.minimize_levels = jb["minimize_levels"].GetBool();
  if (jb.HasMember("fast_mode") && jb["fast_mode"].IsBool())
    ps.fast_mode = jb["fast_mode"].GetBool();
  mockturtle::aig_balance( aig, ps );
}

void apply_rewrite(mockturtle::aig_network& aig, rapidjson::Value const& jb )
{
  mockturtle::rewrite_params ps;

  if (jb.HasMember("cut_limit") && jb["cut_limit"].IsUint())
    ps.cut_enumeration_ps.cut_limit = jb["cut_limit"].GetUint();

  if (jb.HasMember("minimize_truth_table") && jb["minimize_truth_table"].IsBool())
    ps.cut_enumeration_ps.minimize_truth_table = jb["minimize_truth_table"].GetBool();

  if (jb.HasMember("preserve_depth") && jb["preserve_depth"].IsBool())
    ps.preserve_depth = jb["preserve_depth"].GetBool();

  if (jb.HasMember("allow_multiple_structures") && jb["allow_multiple_structures"].IsBool())
    ps.allow_multiple_structures = jb["allow_multiple_structures"].GetBool();

  if (jb.HasMember("allow_zero_gain") && jb["allow_zero_gain"].IsBool())
    ps.allow_zero_gain = jb["allow_zero_gain"].GetBool();

  if (jb.HasMember("use_dont_cares") && jb["use_dont_cares"].IsBool())
    ps.use_dont_cares = jb["use_dont_cares"].GetBool();

  if (jb.HasMember("window_size") && jb["window_size"].IsUint())
    ps.window_size = jb["window_size"].GetUint();

  if (jb.HasMember("verbose") && jb["verbose"].IsBool())
    ps.verbose = jb["verbose"].GetBool();

  mockturtle::xag_npn_resynthesis<mockturtle::aig_network> resyn;
  mockturtle::exact_library_params eps;
  eps.np_classification = false;
  mockturtle::exact_library<mockturtle::aig_network> exact_lib( resyn, eps );

  mockturtle::rewrite(aig, exact_lib, ps);
}

void apply_aig_opto(mockturtle::aig_network& aig, const aig_opto_t& step)
{
  for (int i = 0; i < step.iterations; ++i)
  {
    if (step.algorithm == "balance")
    {
      apply_aig_balancing(aig, *step.params);
    }
    else if (step.algorithm == "rewrite")
    {
      apply_rewrite(aig, *step.params);
    }
    else
    {
      throw std::runtime_error("Unknown algorithm: " + step.algorithm);
    }
  }
}

} // namespace rinox::experiments
