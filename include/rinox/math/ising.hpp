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
#include <assert.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace rinox {

namespace math {

template <typename Data, uint32_t kMaxNumSpins, uint32_t kMaxNumOutputs>
class SpinSystem {
  struct Spin {
    uint32_t name;
    Data info;
    Data care;
    int value;
  };

 public:
  void AddSpin(Data const& info, Data const& care, uint32_t name) {
    int n = spins_.size();
    for (int i = 0; i < n; ++i) {
      int Jp = kitty::count_ones(care & ~info ^ spins_[i].info);
      int Jn = kitty::count_ones(care & info ^ spins_[i].info);
      J[n][i] = (Jp - Jn) / 2;
      J[i][n] = J[n][i];
      std::cout << "J[" << n << "][" << i << "]=" << J[n][i] << std::endl;
    }

    Spin s{name, info, care, 0};
    spins_.push_back(s);
  }

  void clear() { spins_.clear(); }

  int ComputeEnergy() const {
    int res = 0;
    const int n = static_cast<int>(spins_.size());
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        res -= J[i][j] * spins_[i].value * spins_[j].value;
    return res;
  }

  bool GaugeTest() {
    std::cout << "E[0]=" << ComputeEnergy() << std::endl;
    spins_[0].value = 1;
    for (int i = 1; i < spins_.size(); ++i) {
      spins_[i].value = (J[0][i] >= 0) ? 1 : -1;
      std::cout << "Spin " << i << " set to " << spins_[i].value << std::endl;
    }

    std::cout << "E[T]=" << ComputeEnergy() << std::endl;

    for (int i = 0; i < spins_.size() - 1; ++i) {
      for (int j = i + 1; j < spins_.size(); ++j) {
        if (std::signbit(spins_[i].value * spins_[j].value) !=
            std::signbit(J[i][j]))
          return false;
      }
    }

    return true;
  }

  void MinimizeEnergy() {
    bool const trivial = GaugeTest();
    if (trivial) std::cout << "Trivial system detected" << std::endl;
  }

  template <typename Fn>
  void ForeachSpin(Fn&& fn) const {
    for (auto const& spin : spins_) fn(spin);
  }

 private:
  int J[kMaxNumSpins][kMaxNumSpins];
  math::static_vector<Spin, kMaxNumSpins> spins_;
};

}  // namespace math

}  // namespace rinox