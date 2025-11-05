/* rinox: C++ logic network library
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
  \file constexpr_functions.hpp
  \brief Various functions for compile time evaluation

  \author Andrea Costamagna
*/

#pragma once
#include <array>
#include <assert.h>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace rinox
{

namespace math
{

template<typename Data, uint32_t kMaxNumSpins, uint32_t kMaxNumOutputs>
class SpinSystem
{
  struct Spin{
    uint32_t name;
    Data info;
    Data care;
    int value;
  };

public:

  void AddSpin(Data const& info, Data const& care, uint32_t name)
  {
    int n = spins_.size();
    for ( int i = 0; i < n; ++i )
    {
      int Jp = kitty::count_ones( info ^ spins_[i].info );
      int Jn = kitty::count_ones( ~info ^ spins_[i].info );
      J[n][i] = J[i][n] = (Jp-Jn)/2;
    }

    Spin s{name, info, care, 0};
    spins_.push_back(s);
  }

  void clear()
  {
    spins_.clear();
  }

  int ComputeEnergy() const
  {
    int res = 0;
    const int n = static_cast<int>(spins_.size());
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        res -= J[i][j] * spins_[i].value * spins_[j].value;
    return res;
  }

private:
  int J [kMaxNumSpins][kMaxNumSpins];
  math::static_vector<Spin, kMaxNumSpins> spins_;
};

} // namespace math

} // namespace rinox