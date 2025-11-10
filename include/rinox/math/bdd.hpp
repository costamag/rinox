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
  \file bdd.hpp
  \brief Data structure for Binary Decision Trees

  \author Andrea Costamagna
*/

#ifndef RINOX_BDD_HPP_
#define RINOX_BDD_HPP_

#include "math.hpp"

namespace rinox {

namespace math {

template <typename Signal, uint32_t kCapacity>
class BinaryDecisionTree {
 public:
  static constexpr uint32_t kConst0 = 0u;
  static constexpr uint32_t kConst1 = 1u;

 private:
  struct TreeNode {
    Signal v;
    uint32_t t;
    uint32_t e;
  };

 public:
  BinaryDecisionTree() {
    TreeNode const0{0u, 0u, 0u};
    TreeNode const1{1u, 1u, 1u};
    data_.push_back(const0);
    data_.push_back(const1);
  };

  void Reset() { data_.resize(2u); }

  uint32_t AddNode(Signal const& v, uint32_t t, uint32_t e) {
    uint32_t index = data_.size();
    TreeNode n{v, t, e};
    data_.push_back(n);
    return index;
  }

  template <typename Fn>
  void ForeachNode(Fn&& fn) {
    for (auto i = 2u; i < data_.size(); ++i) {
      fn(data_[i], i);
    }
  }

  static_vector<TreeNode, kCapacity + 2u> data_;
};

}  // namespace math

}  // namespace rinox

#endif  // RINOX_BDD_HPP_