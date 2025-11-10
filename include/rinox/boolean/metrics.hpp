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
  \file metrics.h
  \brief Boolean metrics to evaluate Boolean functions

  \author Andrea Costamagna
*/

#ifndef RINOX_BOOLEAN_METRICS_H_
#define RINOX_BOOLEAN_METRICS_H_

#include <kitty/static_truth_table.hpp>
#include <kitty/ternary_truth_table.hpp>

namespace rinox {
namespace boolean {

template <typename TT>
uint32_t GiniImpurity(TT const& bits, TT const& care) {
  uint32_t const num_zeros = kitty::count_ones(care & bits);
  uint32_t const num_cares = kitty::count_ones(care);
  return num_zeros * (num_cares - num_zeros);
}

template <typename TT>
uint32_t GiniImpurity(kitty::ternary_truth_table<TT> const& func) {
  return GiniImpurity(func._bits, func._care);
}

template <typename TT>
double MutualInformation(TT const& bits_x, TT const& bits_y, TT const& care) {
  using u64 = uint64_t;
  const u64 N = kitty::count_ones(care);
  if (N == 0) return 0.0;

  const u64 cx0 = kitty::count_ones(care & ~bits_x);
  const u64 cy0 = kitty::count_ones(care & ~bits_y);

  const u64 c00 = kitty::count_ones(care & (~bits_x) & (~bits_y));
  const u64 c01 = kitty::count_ones(care & (~bits_x) & bits_y);
  const u64 c10 = kitty::count_ones(care & bits_x & (~bits_y));
  const u64 c11 = N - (c00 + c01 + c10);

  auto H = [N](const u64* c, int k) {
    double s = 0.0;
    for (int i = 0; i < k; ++i)
      if (c[i]) s += c[i] * std::log2(static_cast<double>(c[i]));
    return std::log2(static_cast<double>(N)) - s / static_cast<double>(N);
  };

  const u64 mx[2] = {cx0, N - cx0};
  const u64 my[2] = {cy0, N - cy0};
  const u64 mxy[4] = {c00, c01, c10, c11};

  double mi = H(mx, 2) + H(my, 2) - H(mxy, 4);
  if (mi < 0.0 && mi > -1e-12) mi = 0.0;  // clamp tiny negatives
  return mi;
}

template <typename TT>
double Smoothness(TT const& tt) {
  double res = 0;
  size_t const num_vars = tt.num_vars();
  for (size_t i = 0; i < num_vars; ++i) {
    TT const tt0 = kitty::cofactor0(tt, i);
    TT const tt1 = kitty::cofactor1(tt, i);
    res += kitty::count_ones(tt0 ^ ~tt1) / (pow(2, num_vars) * num_vars);
  }
  return res;
}

}  // namespace boolean
}  // namespace rinox

#endif  // RINOX_SYNTHESIS_LUT_SYNTH_METHOD_H_
