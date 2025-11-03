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

#ifndef RINOX_DEPENDENCY_ORACLE_PIVOTS_H_
#define RINOX_DEPENDENCY_ORACLE_PIVOTS_H_

#include <cstdint>

namespace rinox {
namespace dependency {
namespace window {

/*\brief Oracle for greedy selecton of dependency cuts */
template<typename Ntk, typename Manager, typename Simulator, uint32_t kMaxNumInputs>
class OraclePivots{

  static constexpr uint32_t kMaxNumMasks = 1u << kMaxNumInputs;
  static constexpr uint32_t kMaxNumDivs = Manager::divs_capacity;
  static constexpr uint32_t kMaxNumOutputs = Ntk::max_num_outputs;
  using TruthTable = typename Simulator::TruthTable;
  using Spfds = boolean::spfd_manager<TruthTable, kMaxNumMasks>;
  using Signal = typename Ntk::signal;
  using Cut = math::static_vector<Signal, kMaxNumInputs>;

public:
  OraclePivots(Ntk&ntk, Manager&manager, Simulator&simulator):ntk_(ntk),manager_(manager),simulator_(simulator){}
  void Init();
  bool IsCovered();
  uint32_t Evaluate(TruthTable const&);
  void Update(TruthTable const&);

  std::optional<Cut> GetNext(uint32_t);
  
private:
  Ntk & ntk_;
  Manager & manager_;
  Simulator & simulator_;
  Spfds spfds_;
  Cut cut_;
  uint32_t pivot_;
  math::static_vector<Signal, kMaxNumDivs> sorted_divs_;
};

template<typename Ntk, typename Manager, typename Simulator, uint32_t kMaxNumInputs>
void OraclePivots<Ntk, Manager, Simulator, kMaxNumInputs>::Init()
{
  std::vector<TruthTable const*> bits;
  ntk_.foreach_output(manager_.get_pivot(), [&](auto const& f){
    bits.push_back(&simulator_.Get(f));
  });
  spfds_.init( bits);

  cut_.clear();
  pivot_=0;
  sorted_divs_=manager_.get_divs();
  std::stable_sort( sorted_divs_.begin(), sorted_divs_.end(), [&]( auto const& a, auto const& b ){
    return spfds_.evaluate( simulator_.Get(a) )<spfds_.evaluate( simulator_.Get(b) );
  } );
}

template<typename Ntk, typename Manager, typename Simulator, uint32_t kMaxNumInputs>
auto OraclePivots<Ntk, Manager, Simulator, kMaxNumInputs>::GetNext(uint32_t kMaxInps) -> std::optional<Cut>
{
  cut_.clear();
  spfds_.reset();
  // iterate until success
  while(!spfds_.is_covered() && !spfds_.is_saturated() && pivot_ < sorted_divs_.size()){
    Signal pivot = sorted_divs_[pivot_];
    cut_.clear();
    spfds_.reset();
    cut_.push_back( pivot );
    spfds_.update( simulator_.Get(pivot) );
    while(!spfds_.is_covered() && !spfds_.is_saturated() && (cut_.size() < kMaxInps))
    {
      uint32_t best_cost = spfds_.get_num_edges();
      std::optional<Signal> best;
      for ( uint32_t d = pivot_+1; d<sorted_divs_.size(); ++d )
      {
        uint32_t cost = spfds_.evaluate( simulator_.Get(sorted_divs_[d]) );
        if ( cost < best_cost )
        {
          best_cost = cost;
          best = sorted_divs_[d];
        }
      }
      if ( best )
      {
        cut_.push_back(*best);
        spfds_.update( simulator_.Get(*best) );
      }
      else
        break;
    }
    pivot_++;
    if ( spfds_.is_covered() )
    {
      return cut_;
    }
    else
    {
      cut_.clear();
      spfds_.reset();
    }
  }
  return std::nullopt;
}

template<typename Ntk, typename Manager, typename Simulator, uint32_t kMaxNumInputs>
uint32_t OraclePivots<Ntk, Manager, Simulator, kMaxNumInputs>::Evaluate(TruthTable const& tt)
{
  return spfds_.evaluate(tt);
}

template<typename Ntk, typename Manager, typename Simulator, uint32_t kMaxNumInputs>
void OraclePivots<Ntk, Manager, Simulator, kMaxNumInputs>::Update(TruthTable const& tt)
{
  spfds_.update(tt);
}

template<typename Ntk, typename Manager, typename Simulator, uint32_t kMaxNumInputs>
bool OraclePivots<Ntk, Manager, Simulator, kMaxNumInputs>::IsCovered()
{
  return spfds_.is_covered();
}

}  // namespace rinox
}  // namespace dependency
}  // namespace window

#endif  // RINOX_DEPENDENCY_ORACLES_H_