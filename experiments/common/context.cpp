#include "context.hpp"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

namespace rinox::experiments
{

static inline std::string trim( std::string s )
{
  auto is_space = []( unsigned char c ) { return std::isspace( c ) != 0; };
  s.erase( s.begin(), std::find_if( s.begin(), s.end(), [&]( unsigned char c ) { return !is_space( c ); } ) );
  s.erase( std::find_if( s.rbegin(), s.rend(), [&]( unsigned char c ) { return !is_space( c ); } ).base(), s.end() );
  return s;
}

static fs::path default_root()
{
  return fs::path(RINOX_SOURCE_DIR);
}

std::string get_extension(std::string ext)
{
  for (auto& c : ext)
    c = char(std::tolower(unsigned(c)));
  if (ext == "aiger" || ext == "aig")
    return ".aig";
  if (ext == "aag")
    return ".aag";
  if (ext == "verilog" || ext == "v")
    return ".v";
  if (ext == "blif")
    return ".blif";
  if (!ext.empty() && ext.front() != '.')
    ext.insert(ext.begin(), '.');
  return ext;
}

void collect_extension(input_t& input, std::string const& ext, lorina::diagnostic_engine* diag)
{
  input.extension = get_extension(ext);
}

void collect_excluded(input_t& input, const rapidjson::Value& v, lorina::diagnostic_engine* diag)
{
  std::set<std::string> excluded;
  for ( auto it = v.Begin(); it != v.End(); ++it )
    excluded.insert(std::string( it->GetString() ));
  input.excluded = excluded;
}

std::string get_suite(std::string const& name, lorina::diagnostic_engine* diag)
{
  static const std::set<std::string> epfl{
    "adder","bar","div","hyp","log2","max","multiplier",
    "sin","sqrt","square","arbiter","cavlc","ctrl","dec",
    "i2c","int2float","mem_ctrl","priority","router","voter" };

  static const std::set<std::string> iwls{
    "ac97_ctrl","aes_core","des_area","des_perf","DMA",
    "DSP","ethernet","iwls05_i2c","leon2","leon3_opt",
    "leon3","leon3mp","iwls05_mem_ctrl","netcard",
    "pci_bridge32","RISC","sasc","simple_spi","spi",
    "ss_pcm","systemcaes","systemcdes","tv80","usb_funct","usb_phy" };

  static const std::set<std::string> iscas{
    "c17","c432","c499","c880","c1355","c1908",
    "c2670","c3540","c5315","c6288","c7552" };

  if (iscas.count(name)) return "iscas";
  if (iwls.count(name))  return "iwls";
  if (epfl.count(name))  return "epfl";

  rinox::diagnostics::REPORT_DIAG(diag, lorina::diagnostic_level::fatal, "`{}` benchmark not found", name);
  return "";
}

bool is_benchmark_suite(std::string const& name)
{
  return name == "epfl" || name == "iscas" || name == "iwls";
}

fs::path default_suite_list_file(const fs::path& suite_dir, const std::string& suite_name)
{
  fs::path p = suite_dir / (suite_name + ".suite");
  if (std::filesystem::exists(p)) return p;
  throw std::runtime_error("Suite list file not found: " + p.string());
}

bool read_text_list_file(const fs::path& path, std::vector<std::string>& out_lines)
{
  std::ifstream in(path);
  if (!in) return false;

  std::string line;
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;
    out_lines.push_back(line);
  }
  return true;
}

void collect_benchmark_named(input_t& input, std::string const& name, lorina::diagnostic_engine* diag)
{
  std::string suite = get_suite(name, diag);
  if (suite.empty()) return;

  std::string subdir;
  if (input.extension == ".aig")
    subdir = "aiger";
  else
  {
    rinox::diagnostics::REPORT_DIAG(diag, lorina::diagnostic_level::fatal, "Input extension `{}` not supported", input.extension.c_str());
    return;
  }

  auto const file_name = name + input.extension;
  fs::path path = default_root() / "benchmarks" / suite / subdir / file_name;
  if (!std::filesystem::exists(path)) {
    rinox::diagnostics::REPORT_DIAG(diag, lorina::diagnostic_level::fatal, "`{}` benchmark not found", path.c_str());
    return;
  }

  input.bench_names.push_back(name);
  input.bench_paths.push_back(path);
}

void collect_benchmark_suite(input_t& input, std::string const& suite, lorina::diagnostic_engine* diag)
{
  fs::path suite_dir = default_root() / "benchmarks" / suite;
  fs::path list_path = default_suite_list_file(suite_dir, suite);

  std::string subdir;
  if (input.extension == ".aig")
    subdir = "aiger";
  else {
    rinox::diagnostics::REPORT_DIAG(diag, lorina::diagnostic_level::fatal, "Input extension `{}` not supported", input.extension.c_str());
    return;
  }

  std::vector<std::string> listed;
  if (!read_text_list_file(list_path, listed)) {
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::warning, "Could not read suite list file");
    return;
  }

  for (auto const& base : listed) {
    if (!input.excluded || ((*input.excluded).count(base) == 0)) {
      std::string const file_name = base + input.extension;
      fs::path path = suite_dir / subdir / file_name;
      if (std::filesystem::exists(path))
      {
        input.bench_names.push_back(base);
        input.bench_paths.push_back(path);
      }
    }
  }
}

void collect_benchmark(input_t& input, std::string const& name, lorina::diagnostic_engine* diag)
{
  if (is_benchmark_suite(name))
    collect_benchmark_suite(input, name, diag);
  else
    collect_benchmark_named(input, name, diag);
}

void collect_benchmarks(input_t& input, const rapidjson::Value& arr, lorina::diagnostic_engine* diag)
{
  for (auto const& val : arr.GetArray())
    if (val.IsString())
      collect_benchmark(input, val.GetString(), diag);
}

input_t parse_context_input(const rapidjson::Value& jb, lorina::diagnostic_engine* diag)
{
  input_t input;

  if (jb.HasMember("extension") && jb["extension"].IsString())
    collect_extension(input, jb["extension"].GetString(), diag );
  else {
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::fatal, "Required field `extension` is missing");
    return input;
  }

  if (jb.HasMember("excluded") && jb["excluded"].IsArray())
    collect_excluded(input, jb["excluded"], diag);

  if (jb.HasMember("benchmarks") && jb["benchmarks"].IsArray())
    collect_benchmarks(input, jb["benchmarks"], diag);
  else {
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::fatal, "Required field `benchmarks` is missing");
    return input;
  }

  return input;
}

technology_t parse_context_technology( const rapidjson::Value& jb, lorina::diagnostic_engine* diag )
{
  technology_t tech;
  std::string file_type;
  if ( jb.HasMember("extension") && jb["extension"].IsString() )
  {
    if ( (jb["extension"] == "genlib")||(jb["extension"] == ".genlib") )
    {
      tech.extension = ".genlib";
      file_type = "genlib";
    }
    else
    {
      rinox::diagnostics::REPORT_DIAG(diag, lorina::diagnostic_level::fatal, "`{}` is not a supported technology format", jb["extension"].GetString());
      return tech;
    }
  }
  else
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::note, "Extension not specified. Assumed `.genlib`");


  if ( jb.HasMember("name") && jb["name"].IsString() )
  {
    std::string const name = jb["name"].GetString();
    std::string const file = name + tech.extension;
    tech.path = default_root() / "techlib" / file_type / file;
    if ( !std::filesystem::exists( tech.path ) )
    {
      rinox::diagnostics::REPORT_DIAG(diag, lorina::diagnostic_level::fatal, "`{}` library not found", tech.path);
    }
  }
  else
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::note, "Extension not specified. Assumed `.genlib`");
  return tech;
}

output_t parse_context_output( const rapidjson::Value& jb, lorina::diagnostic_engine* diag )
{
  output_t out;
  std::string file_type;
  if ( jb.HasMember("path") && jb["path"].IsString() )
  {
    if ( !std::filesystem::exists( jb["path"].GetString() ) )
    {
      rinox::diagnostics::REPORT_DIAG(diag, lorina::diagnostic_level::fatal, "`{}` directory not found", jb["path"].GetString());
      return out;
    }
    else
      out.path = jb["path"].GetString();
  }

  return out;
}

context_t load_context(const rapidjson::Document& doc, lorina::diagnostic_engine* diag)
{
  context_t context;

  if (doc.HasMember("context") && doc["context"].IsObject()) {
    const auto& jb = doc["context"];
    
    // technology section
    if (jb.HasMember("technology") && jb["technology"].IsObject())
      context.technology = parse_context_technology( jb["technology"], diag );
    else {
      rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::fatal, "A valid technology must be specified");
      return context;
    }

    // input section
    if (jb.HasMember("input") && jb["input"].IsObject())
      context.input = parse_context_input(jb["input"], diag);

    // output section
    if (jb.HasMember("output") && jb["output"].IsObject())
      context.output = parse_context_output(jb["output"], diag);

  } else {
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::fatal, ".config files must contain the context field");
    return context;
  }

  if (!is_specified(context)) {
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::fatal, "Context not specified after parsing");
    return context;
  }

  return context;
}

std::optional<std::vector<mockturtle::gate>> load_gates( context_t const& ctx, lorina::diagnostic_engine * diag )
{
  std::vector<mockturtle::gate> gates;
  std::ifstream in( ctx.technology.path );

  if ( lorina::read_genlib( in, mockturtle::genlib_reader( gates ) ) != lorina::return_code::success )
  {
    rinox::diagnostics::REPORT_DIAG_RAW(diag, lorina::diagnostic_level::fatal, "Failed loading library");
    return std::nullopt;
  }
  return gates;
}

} // namespace rinox::experiments
