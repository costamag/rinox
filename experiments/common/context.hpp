#pragma once
#include <rinox/traits.hpp>
#include <rinox/diagnostics/diagnostics.hpp>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <mockturtle/io/genlib_reader.hpp>
#include <filesystem>
#include <set>
#include <string>
#include <vector>

namespace rinox::experiments
{

namespace fs = std::filesystem;

template<typename T>
std::enable_if_t<std::is_fundamental_v<T>, bool>
is_specified(const T& val) {
  return val != T{};
}

template<typename T>
std::enable_if_t<rinox::traits::is_container_v<T>, bool>
is_specified(const T& val) {
  return !val.empty();
}

template<typename T>
std::enable_if_t<rinox::traits::is_optional_v<T>, bool>
is_specified(const T&) {
  return true; // presence is optional
}

template<typename T>
std::enable_if_t<
  !std::is_fundamental_v<T> &&
  !rinox::traits::is_container_v<T> &&
  !rinox::traits::is_optional_v<T>, bool>
is_specified(const T& obj) {
  bool result = true;
  inspect_fields(obj, [&](const char*, const auto& field) {
    if constexpr (!rinox::traits::is_optional_v<std::decay_t<decltype(field)>>) {
      result = result && is_specified(field);
    }
  });
  return result;
}

struct input_t {
  std::string extension;
  std::vector<std::string> bench_paths;
  std::vector<std::string> bench_names;
  std::optional<std::set<std::string>> excluded;
};

template<typename Visitor>
void inspect_fields(const input_t& input, Visitor&& visit) {
  visit("extension", input.extension);
  visit("bench_paths", input.bench_paths);
  visit("bench_names", input.bench_names);
}

struct output_t {
  std::string path;
};

template<typename Visitor>
void inspect_fields(const output_t& output, Visitor&& visit) {
  visit("path", output.path);
}

struct technology_t {
  std::string extension;
  std::string path;
};

template<typename Visitor>
void inspect_fields(const technology_t& tech, Visitor&& visit) {
  visit("extension", tech.extension);
  visit("path", tech.path);
}


struct context_t {
  technology_t technology;
  std::optional<input_t> input;
  std::optional<output_t> output;
};

template<typename Visitor>
void inspect_fields(const context_t& ctx, Visitor&& visit) {
  visit("technology", ctx.technology); // optional
  visit("input", ctx.input);
  visit("output", ctx.output);         // optional
}

context_t load_context( const rapidjson::Document& doc, lorina::diagnostic_engine* diag = nullptr );

template<typename Fn>
bool foreach_benchmark( context_t const& ctx, Fn && fn )
{
  auto const& bench_paths = (*ctx.input).bench_paths;
  auto const& bench_names = (*ctx.input).bench_names;

  for ( auto i = 0u; i < bench_names.size(); ++i )
  {
    if ( !fn( bench_paths[i], bench_names[i] ) )
      return false;
  }
  return true;
}

std::optional<std::vector<mockturtle::gate>> load_gates( context_t const& ctx, lorina::diagnostic_engine * diag );

} // namespace rinox::experiments
