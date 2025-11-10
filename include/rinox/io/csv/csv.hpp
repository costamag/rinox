/* rinox: C++ parsing library
 * Copyright (C) 2025 EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file json_parser.hpp
  \brief Implements simplistic Json parser

  \author Andrea Costamagna
*/

#ifndef RINOX_IO_CSV_HPP_
#define RINOX_IO_CSV_HPP_

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <locale>
#include <map>
#include <string>
#include <vector>

namespace rinox {

namespace io {

namespace csv {

// Writes a CSV with columns given by name -> vector<T>.
// If vectors have different lengths, missing cells are left empty.
// Good for numeric types; doubles are written with high precision.
template <typename T>
void WriteCsv(const std::string& filepath,
              const std::map<std::string, std::vector<T>>& columns,
              bool include_index = true, const std::string& index_name = "i") {
  // Compute row count as max column size
  std::size_t rows = 0;
  for (const auto& [_, col] : columns) rows = std::max(rows, col.size());

  std::ofstream out(filepath, std::ios::binary);
  if (!out) throw std::runtime_error("write_csv: cannot open " + filepath);

  // Ensure dot as decimal separator regardless of system locale
  out.imbue(std::locale::classic());

  // High precision for floating point
  out << std::setprecision(17) << std::boolalpha;

  auto write_cell = [&](const std::string& s, bool last) {
    // Numbers don’t need escaping; names might. Quote if needed.
    bool needs_quotes = s.find_first_of(",\"\n\r") != std::string::npos;
    if (needs_quotes) {
      out << '"';
      for (char c : s) out << (c == '"' ? "\"\"" : std::string(1, c));
      out << '"';
    } else {
      out << s;
    }
    if (!last) out << ',';
  };

  // Header
  bool first = true;
  if (include_index) {
    write_cell(index_name, false);
    first = false;
  }
  std::size_t ci = 0, ccount = columns.size();
  for (const auto& [name, _] : columns) {
    write_cell(name, ci + 1 == ccount);
    ++ci;
  }
  out << '\n';

  // Rows
  for (std::size_t r = 0; r < rows; ++r) {
    if (include_index) write_cell(std::to_string(r), columns.empty());
    ci = 0;
    for (const auto& [_, col] : columns) {
      if (r < col.size()) {
        std::ostringstream ss;
        ss.imbue(std::locale::classic());
        ss << std::setprecision(17) << col[r];
        write_cell(ss.str(), ci + 1 == ccount);
      } else {
        write_cell("", ci + 1 == ccount);
      }
      ++ci;
    }
    out << '\n';
  }
}

// Convenience overload for doubles with vector names inferred from order
inline void WriteCsv(
    const std::string& filepath,
    const std::vector<std::pair<std::string, std::vector<double>>>& named_cols,
    bool include_index = true, const std::string& index_name = "i") {
  std::map<std::string, std::vector<double>> m(named_cols.begin(),
                                               named_cols.end());
  WriteCsv(filepath, m, include_index, index_name);
}

}  // namespace csv
}  // namespace io
}  // namespace rinox
#endif  // RINOX_IO_CSV_HPP_