/* rinox: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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

/* rinox: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
 * ...license...
 */

/*!
  \file lut_synth_method.h
  \brief Synthesis method enum for LUT flows
  \author Andrea Costamagna
*/

#ifndef RINOX_SYNTHESIS_LUT_SYNTH_METHOD_H_
#define RINOX_SYNTHESIS_LUT_SYNTH_METHOD_H_

#include "shannon.hpp"
#include "spin_dec.hpp"
#include <cstdint>  // uint8_t

namespace rinox {
namespace synthesis {
namespace lut {

enum class LutSynthMethod : uint8_t {
  kShannonSolo,
  kSpinDecSolo
};

template<LutSynthMethod M, class Ntk, class Manager, class Simulator, class Params>
struct LutSynthDispatch;  // primary template (no definition)

// Shannon solo
template<class Ntk, class Manager, class Simulator, class Params>
struct LutSynthDispatch<LutSynthMethod::kShannonSolo, Ntk, Manager, Simulator, Params> {
  using type = ShannonDecSolo<Ntk, Manager, Simulator, Params>;
};

// Ising solo
template<class Ntk, class Manager, class Simulator, class Params>
struct LutSynthDispatch<LutSynthMethod::kSpinDecSolo, Ntk, Manager, Simulator, Params> {
  using type = SpinDecSolo<Ntk, Manager, Simulator, Params>;
};

}  // namespace lut
}  // namespace synthesis
}  // namespace rinox

#endif  // RINOX_SYNTHESIS_LUT_SYNTH_METHOD_H_
