/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
  \file verilog_reader.hpp
  \brief Lorina reader for VERILOG files

  \author Andrea Costamagna
  \author Heinz Riener
  \author Marcel Walter
  \author Mathias Soeken
  \author Siang-Yun (Sonia) Lee
*/

#pragma once

#include <rinox/diagnostics/diagnostics.hpp>
#include "../utils/reader.hpp"
#include "verilog_parser.hpp"
#include "write_verilog.hpp"

namespace rinox
{

namespace io
{

namespace verilog
{

#ifndef DIAG_HERE_FILE
#define DIAG_HERE_FILE __FILE__
#endif

#ifndef DIAG_HERE_LINE
#define DIAG_HERE_LINE __LINE__
#endif

/*! \brief Reader function for VERILOG format.
 *
 * Reads a simplistic VERILOG format from a stream and invokes a callback
 * method for each parsed primitive and each detected parse error.
 *
 * \param in Input stream
 * \param reader A VERILOG reader with callback methods invoked for parsed primitives
 * \param diag An optional diagnostic engine with callback methods for parse errors
 * \return Success if parsing has been successful, or parse error if parsing has failed
 */
template<typename Ntk>
[[nodiscard]] inline lorina::return_code read_verilog( std::istream& in, const reader<Ntk>& reader, lorina::diagnostic_engine* diag = nullptr )
{
  verilog_parser parser( in, reader, diag );
  auto result = parser.parse_modules();
  if ( !result )
  {
    return lorina::return_code::parse_error;
  }
  else
  {
    return lorina::return_code::success;
  }
}

/*! \brief Reader function for VERILOG format.
 *
 * Reads a simplistic VERILOG format from a file and invokes a callback
 * method for each parsed primitive and each detected parse error.
 *
 * \param filename Name of the file
 * \param reader A VERILOG reader with callback methods invoked for parsed primitives
 * \param diag An optional diagnostic engine with callback methods for parse errors
 * \return Success if parsing has been successful, or parse error if parsing has failed
 */
template<typename Ntk>
[[nodiscard]] inline lorina::return_code read_verilog( const std::string& filename, const reader<Ntk>& reader, lorina::diagnostic_engine* diag = nullptr )
{
  std::ifstream in( lorina::detail::word_exp_filename( filename ), std::ifstream::in );
  if ( !in.is_open() )
  {
    rinox::diagnostics::REPORT_DIAG( diag, lorina::diagnostic_level::fatal, "failed to open file `{}`", filename.c_str() );
    return lorina::return_code::parse_error;
  }
  else
  {
    auto const ret = read_verilog( in, reader, diag );
    in.close();
    if ( ret != lorina::return_code::success )
      rinox::diagnostics::REPORT_DIAG( diag, lorina::diagnostic_level::fatal, "failed to read the verilog file `{}`", filename.c_str() );

    return ret;
  }
}

} /* namespace verilog */

} /* namespace io */

} // namespace rinox