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
  \file struct.hpp
  \brief Structural resinthesis of LUT networks

  \author Andrea Costamagna
*/

#ifndef RINOX_WINDOWING_WINDOW_TAGS_HPP_
#define RINOX_WINDOWING_WINDOW_TAGS_HPP_

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <mockturtle/utils/node_map.hpp>
#include <optional>
#include <queue>
#include <stack>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../math/math.hpp"
#include "../network/visualize.hpp"

namespace rinox {

namespace windowing {

struct WindowTags {
  enum class tag : uint32_t {
    INPS = 1,
    DIVS = 2,
    MFFC = 4,
    TFOS = 8,
    MASK = 0xF
  };

  static constexpr uint32_t num_bits_tag = 4u;  // low bits
  static constexpr uint32_t role_mask = static_cast<uint32_t>(tag::MASK);
  static constexpr uint32_t id_mask = ~role_mask;  // high bits

  uint32_t value{0};

  WindowTags() = default;
  explicit WindowTags(tag t) : value(static_cast<uint32_t>(t)) {}

  static WindowTags FromId(uint32_t id) {
    if (id > (UINT32_MAX >> num_bits_tag))
      throw std::out_of_range("color_t id too large");
    return WindowTags((id << num_bits_tag), true);
  }

  operator uint32_t() const { return value; }

  void Add(tag t) { value = (value & id_mask) | static_cast<uint32_t>(t); }
  void Remove(tag) { value &= id_mask; }  // remove role bits entirely

  bool In(tag t) const {
    uint32_t tag_bits = value & role_mask;
    if (tag_bits == 0) return false;
    if (tag_bits == static_cast<uint32_t>(t)) return true;
    if (tag_bits & static_cast<uint32_t>(t)) {
      throw std::logic_error("Node has incompatible colors");
    }
    return false;
  }

  uint32_t VisitedId() const { return value >> num_bits_tag; }
  void SetVisited(uint32_t id) {
    value = (id << num_bits_tag) | (value & role_mask);
  }

 private:
  explicit WindowTags(uint32_t raw, bool) : value(raw) {}
};

}  // namespace windowing
}  // namespace rinox

#endif