/* lorina: C++ parsing library
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
  \file verilog.hpp
  \brief Implements simplistic Verilog parser

  \author Andrea Costamagna
*/

#pragma once

#include <fmt/color.h>
#include <fmt/format.h>
#include <lorina/common.hpp>
#include <lorina/diagnostics.hpp>
#include <string>

namespace rinox
{

namespace diagnostics
{

#ifndef DIAG_HERE_FILE
#define DIAG_HERE_FILE __FILE__
#endif

#ifndef DIAG_HERE_LINE
#define DIAG_HERE_LINE __LINE__
#endif

template<typename... Args>
inline void report_diagnostic( lorina::diagnostic_engine* diag,
                               lorina::diagnostic_level level,
                               const char* fmt,
                               const char* file, int line,
                               Args&&... args )
{
  if ( !diag )
    return;

  // primary message
  auto id = diag->create_id( level, fmt );
  auto rep = diag->report( id );

  // call rep.add_argument(arg) for each arg in the pack
  (void)std::initializer_list<int>{
      ( (void)rep.add_argument( std::forward<Args>( args ) ), 0 )... };

  // clickable location as a separate note
  if ( level != lorina::diagnostic_level::note )
  {
    auto note_code = diag->create_id( level, "  ↪ {}:{}" );
    diag->report( note_code ).add_argument( std::string( file ) ).add_argument( std::to_string( line ) );
  }
}

template<typename... Args>
inline void report_diagnostic_raw( 
                                   lorina::diagnostic_engine* diag,
                                   lorina::diagnostic_level level,
                                   const char* fmt,
                                   Args&&... args )
{
  if ( !diag )
    return;

  // primary message
  auto id = diag->create_id( level, fmt );
  auto rep = diag->report( id );

  // call rep.add_argument(arg) for each arg in the pack
  (void)std::initializer_list<int>{
      ( (void)rep.add_argument( std::forward<Args>( args ) ), 0 )... };
}

#define REPORT_DIAG( diag, level, fmt, ... ) \
  report_diagnostic( ( diag ), ( level ), ( fmt ), __FILE__, __LINE__, ##__VA_ARGS__ )

#define REPORT_DIAG_RAW( diag, level, fmt, ... ) \
  report_diagnostic_raw( ( diag ), ( level ), ( fmt ), ##__VA_ARGS__ )

/*! \brief A consumer for diagnostics. */
class text_diagnostics : public lorina::diagnostic_consumer
{
public:
  text_diagnostics() = default;
  virtual ~text_diagnostics() = default;

  /*! \brief Handle diagnostic.
   *
   * \param level Severity level
   * \param message Diagnostic message
   */
  void handle_diagnostic( lorina::diagnostic_level level, std::string const& message ) const override
  {
    switch ( level )
    {
    case lorina::diagnostic_level::ignore:
      break;
    case lorina::diagnostic_level::note:
    {
      fmt::print( stdout, fmt::emphasis::bold | fg( fmt::terminal_color::bright_green ), "[i]" );
      fmt::print( stdout, fmt::emphasis::bold | fg( fmt::terminal_color::bright_green ), " {}\n", message );
    }
    break;
    case lorina::diagnostic_level::remark:
    {
      fmt::print( stderr, fmt::emphasis::bold | fg( fmt::terminal_color::bright_green ), "[I]" );
      fmt::print( stderr, fmt::emphasis::bold | fg( fmt::terminal_color::bright_green ), " {}\n", message );
    }
    break;
    case lorina::diagnostic_level::warning:
    {
      fmt::print( stderr, fmt::emphasis::bold | fg( fmt::terminal_color::bright_magenta ), "[w]" );
      fmt::print( stderr, fmt::emphasis::bold | fg( fmt::terminal_color::bright_magenta ), " {}\n", message );
    }
    break;
    case lorina::diagnostic_level::error:
    {
      fmt::print( stderr, fmt::emphasis::bold | fg( fmt::terminal_color::bright_red ), "[e]" );
      fmt::print( stderr, fmt::emphasis::bold | fg( fmt::terminal_color::bright_red ), " {}\n", message );
    }
    break;
    case lorina::diagnostic_level::fatal:
    default:
    {
      fmt::print( stderr, fmt::emphasis::bold | fg( fmt::terminal_color::bright_red ), "[E]" );
      fmt::print( stderr, fmt::emphasis::bold | fg( fmt::terminal_color::bright_red ), " {}\n", message );
    }
    break;
    }
  }
};

} // namespace diagnostics

} // namespace rinox