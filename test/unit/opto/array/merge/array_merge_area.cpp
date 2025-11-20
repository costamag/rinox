#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <kitty/kitty.hpp>
#include <kitty/static_truth_table.hpp>

#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <rinox/network/network.hpp>
#include <rinox/opto/array/merge/array_merge_area.hpp>


TEST_CASE( "Merging a 2-input node with a 1-input node", "[merge_array]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 2>;
  using signal = typename Ntk::signal;

  Ntk ntk;
  std::vector<signal> pis;
  for (int i=0; i<4; ++i)
    pis.push_back(ntk.create_pi());

  kitty::dynamic_truth_table tt2(2u);
  kitty::dynamic_truth_table tt1(1u);
  kitty::create_from_binary_string(tt2, "1001");
  kitty::create_from_binary_string(tt1, "01");
  auto const f1 = ntk.create_node({pis[0], pis[1]}, tt2);
  auto const f2 = ntk.create_node({pis[3]}, tt1);
  ntk.create_po(f1);
  ntk.create_po(f2);
  using Params = rinox::opto::array::MergeAreaParamsDefault;
  using Stats = rinox::opto::array::MergeAreaStatsDefault;
  Params ps;
  Stats st;
  CHECK(ntk.area()==2);
  rinox::opto::array::MergeArea<Ntk, Params, Stats>(ntk, ps, &st);
  CHECK(ntk.area()==1);
  st.Report();
  
}



TEST_CASE( "Merging a 2-input node with a 1-input node in cascade", "[merge_array]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, 2>;
  using signal = typename Ntk::signal;

  Ntk ntk;
  std::vector<signal> pis;
  for (int i=0; i<4; ++i)
    pis.push_back(ntk.create_pi());

  kitty::dynamic_truth_table tt2(2u);
  kitty::dynamic_truth_table tt1(1u);
  kitty::create_from_binary_string(tt2, "1001");
  kitty::create_from_binary_string(tt1, "01");
  auto const f2 = ntk.create_node({pis[3]}, tt1);
  auto const f1 = ntk.create_node({pis[0], f2}, tt2);
  ntk.create_po(f1);
  ntk.create_po(f2);
  using Params = rinox::opto::array::MergeAreaParamsDefault;
  using Stats = rinox::opto::array::MergeAreaStatsDefault;
  Params ps;
  Stats st;
  CHECK(ntk.area()==2);
  rinox::opto::array::MergeArea<Ntk, Params, Stats>(ntk, ps, &st);
  CHECK(ntk.area()==1);
  st.Report();
  
}
