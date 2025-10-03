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
  \file equivalence_checking.hpp
  \brief Combinational equivalence checking

  \author Andrea Costamagna
*/

#ifndef RINOX_EQUIVALENCE_CHECKING_HPP_
#define RINOX_EQUIVALENCE_CHECKING_HPP_

#include <cstdint>
#include <iostream>
#include <vector>

#include "cleanup.hpp"
#include "functional_reduction.hpp"
#include "../traits.hpp"
#include "../utils/include/percy.hpp"
#include "../utils/stopwatch.hpp"
#include "../networks/klut.hpp"
#include "cnf.hpp"

#include <bill/sat/interface/common.hpp>
#include <bill/sat/interface/abc_bsat2.hpp>

#include <fmt/format.h>

namespace rinox{
namespace windowing {


template<uint32_t kNumInputs>
class CounterExamples{
public:
  void SetBit(uint32_t i){ kitty::set_bit( bits_[i], count_ ); }
  void ClearBit(uint32_t i){ kitty::clear_bit( bits_[i], count_ ); }
  void IncrementCounter(uint32_t i){ count_++; }
  bool IsFull(){ return count_ == 62u; }
  void Reset(){ count_=0; }

private:
  math::static_vector<kitty::static_truth_table<6u>, kNumInputs> bits_;
  uint32_t count_ = 0u;
};

struct CECResynParams
{
  uint32_t random_seed{ 123 };
  uint32_t conflict_limit{ 1000000 };
  bool verbose{ false };
};

struct CECResynStats
{
  std::vector<bool> cex;

  void Report() const
  {
  }
};

template<typename Ntk, typename WindowManager>
class CECResyn{
  static constexpr uint32_t kMaxNumLeaves = WindowManager::kMaxNumLeaves;
  using Params = CECResynParams;
  using Stats = CECResynStats;
  using Signal = typename Ntk::signal;
  using Solver = bill::solver<bill::solvers::bsat2>;
  using CExs = CounterExamples<kMaxNumLeaves>;
  using AddClauseFnType = std::function<void( std::vector<bill::lit_type> const& )>;

  CECResyn(Ntk&, WindowManager&, Params&, Stats&);
  void Restart();
  void AddClause(Signal);
  std::optional<bool> Solve(std::vector<bill::lit_type>);
  std::optional<bool> Solve(std::vector<bill::lit_type>);

private:
  Ntk & ntk_;
  WindowManager & manager_;
  math::incomplete_signal_map<bill::lit_type, Ntk> literals_;
  std::vector<bill::lit_type> input_lits_;
  AddClauseFnType AddClauseFn = [&]( auto const& clause ) { solver_.add_clause( clause ); };

  Solver solver_;
  CExs cexs_;
  Params & ps_;
  Stats & st_;
};

template<typename Ntk, typename WindowManager>
CECResyn<Ntk, WindowManager>::CECResyn(Ntk&ntk, WindowManager & manager, Params&ps, Stats&st)
  : ntk_(ntk), manager_(manager), ps_(ps), st_(st)
{
  literals_[ntk_.get_constant( false )] = bill::lit_type( 0, bill::lit_type::polarities::positive );
  literals_[ntk.get_constant( true )] = bill::lit_type( 0, bill::lit_type::polarities::negative );
  for ( uint32_t lit = 0; lit < kMaxNumLeaves; ++lit )
    input_lits_.push_back( bill::lit_type( i + 1, bill::lit_type::polarities::positive ) );

  Restart();
}

template<typename Ntk, typename WindowManager>
void CECResyn<Ntk, WindowManager>::Restart()
{
  num_invoke_ = 0u;
  solver_.restart();
  solver_.set_random_phase( ps_.random_seed );
  solver_.add_variables( kMaxNumLeaves + 1 );
  solver_.add_clause( { ~literals_[ntk_.get_constant( false )] } );

  literals_.resize();
  manager_.foreach_input([&]( auto const& f, auto i ){
    literals_[f] = input_lits_[i];
  });
  manager_.foreach_divisor([&]( auto const& f, auto i ){
    if ( !manager_.is_input( ntk_.get_node(f) ) ) AddClause(f);
  });
  manager_.foreach_mffc([&]( auto const& n, auto i ){
    ntk_.foreach_output( n, [&]( auto const& f ){ AddClause(f); } );
  });
  manager_.foreach_tfo([&]( auto const& n, auto i ){
    ntk_.foreach_output( n, [&]( auto const& f ){ AddClause(f); } );
  });
}

template<typename Ntk, typename WindowManager>
void CECResyn<Ntk, WindowManager>::AddClause(Signal f)
{
  auto lit = bill::lit_type( solver_.add_variable(), bill::lit_type::polarities::positive );
  std::vector<bill::lit_type> const& child_lits;
  ntk_.foreach_fanin( f, [&](auto const& fi, auto ii){
    child_lits.push_back( literals_[fi] );
  } );
  auto const& function = ntk_.signal_function( f );
  mockturtle::detail::on_function<AddClauseFnType>( lit, child_lits, function, AddClauseFn );
  literals_[f]=lit;
}

template<typename Ntk, typename WindowManager>
std::optional<bool> CECResyn<Ntk, WindowManager>::Solve( std::vector<bill::lit_type> assumptions )
{
  ++num_invoke_;
  auto const res = solver_.solve( assumptions, ps_.conflict_limit );

  if ( res == bill::result::states::satisfiable )
  {
    auto model = solver.get_model().model();
    for ( auto i = 0u; i < kMaxNumLeaves; ++i )
    {
      if ( model.at( i + 1 ) == bill::lbool_type::true_ )
        cexs_.SetBit( i );
      else
        cexs_.ClearBit( i );
    }
    cexs_.IncrementCounter();
    return false;
  }
  else if ( res == bill::result::states::unsatisfiable )
  {
    return true;
  }
  else
  {
    return std::nullopt; /* timeout or something wrong */
  }
}

} // namespace windowing
} // namespace rinox

#endif