#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <kitty/kitty.hpp>

#include <rinox/network/network.hpp>
#include <rinox/synthesis/lut/spin_dec.hpp>
#include <rinox/math/math.hpp>

struct LutSynthParams1 : rinox::synthesis::lut::DefaultLutSynthParams {
  static constexpr uint32_t kNumVarsCube = 6u;
  static constexpr uint32_t kMaxCutsSize = 5u;
  static constexpr uint32_t kMaxInpsSolo = 4u;
  static constexpr uint32_t kMaxInpsDual = 2u;
  static constexpr uint32_t kMaxNumOutputs = 1u;
  static constexpr uint32_t kMaxNumGates = 32u;
  static constexpr bool kExactSuppMin = false;
};

struct ManagerParams : rinox::windowing::default_window_manager_params {
  uint32_t max_num_divisors = 128u;
  uint32_t max_reconvergence = 10u;
  int32_t odc_levels = 0u;
  static constexpr uint32_t max_num_leaves = 6u;
};

template<typename Ntk>
class ToyManager
{
  using Signal = typename Ntk::signal;
  using Node = typename Ntk::node;
public:
  ToyManager(std::vector<Signal> const& divs):
    divs_(divs)
  {}

  private:
    std::vector<Signal> divs_;
};

template<typename Ntk, uint32_t MaxLeaves = 6u>
class ToySimulator
{
  using Signature = kitty::static_truth_table<MaxLeaves>;
  using Signal = typename Ntk::signal;
  using Node = typename Ntk::node;
public:
  ToySimulator(std::vector<Signal> const& divs, std::vector<Signature> const& sims):
    divs_(divs),
    sims_(sims)
  {
    for ( int i = 0; i < divs.size(); ++i )
    {
      uint64_t key = divs_[i].data;
      cache_[key] = i;
    }
  }

    Signature const& Get(Signal s) const
    {
        uint64_t key = s.data;
        auto it = cache_.find(key);
        if ( it != cache_.end() )
        {
        return sims_[ it->second ];
        }
        throw std::runtime_error("Signal not found in simulator");
    }

  private:
    std::vector<Signal> divs_;
    std::vector<Signature> sims_;
    std::unordered_map<uint64_t, int> cache_;
};

template<typename Ntk, typename Signature>
void simulate( Ntk & ntk, typename Ntk::signal f, std::unordered_map<uint64_t, Signature>& map) {
  if (ntk.is_pi(ntk.get_node(f)))
    return;
  ntk.foreach_fanin(f.index, [&](auto const& fi, auto){
    if ( map.find(fi.data) == map.end() )
    {
      simulate(ntk, fi, map);
    }
  });
  std::vector<Signature const *> inps;
  ntk.foreach_fanin(f.index, [&](auto const& fi, auto){
    inps.push_back( &map[fi.data] );
  });
  map[f.data] = ntk.compute( f, inps );
}

TEST_CASE( "Spin-based decomposition - 1 LUT", "[lut_spin_dec]" )
{
  using Ntk = rinox::network::bound_network<rinox::network::design_type_t::ARRAY_BASED, LutSynthParams1::kMaxNumOutputs>;
  using Manager = ToyManager<Ntk>;
  using Simulator = ToySimulator<Ntk>;
  using Signature = kitty::static_truth_table<LutSynthParams1::kNumVarsCube>;
  using Signal = typename Ntk::signal;
  using Node = typename Ntk::node;
  std::unordered_map<uint64_t, Signature> map;

  Ntk ntk;
  std::vector<Signal> divs;
  std::vector<Signature> sims;
  for (int i = 0; i < LutSynthParams1::kNumVarsCube; ++i)
  {
    auto pi = ntk.create_pi();
    divs.push_back( pi );
    sims.emplace_back();
    kitty::create_nth_var( sims[i], i );
    map[pi.data] = sims.back();
  }

  Manager manager(divs);
  Simulator simulator(divs, sims);

  rinox::synthesis::lut::LutSynthStats st;
  rinox::synthesis::lut::SpinDecSolo<Ntk, Manager, Simulator, LutSynthParams1> synth( ntk, manager, simulator, st );

  rinox::math::static_vector<Signal, LutSynthParams1::kMaxCutsSize> cut;
  std::array<Signature, LutSynthParams1::kMaxCutsSize> xs;
  for ( int i = 0; i < LutSynthParams1::kMaxCutsSize; ++i )
  {
    cut.push_back( divs[i] );
    kitty::create_nth_var( xs[i], i );
  }
  Signature func, care;

  std::function<void(Signature &, Signature &)>Check  = [&](Signature & func, Signature & care){
    for ( int i = LutSynthParams1::kMaxCutsSize; i < LutSynthParams1::kNumVarsCube; ++i )
    {
      func = func & care;
      auto func0 = kitty::cofactor0(func, i);
      auto func1 = kitty::cofactor1(func, i);
      auto care0 = kitty::cofactor0(care, i);
      auto care1 = kitty::cofactor1(care, i);
      func = func0 | func1;
      care = care0 | care1;
    }
    rinox::math::static_vector<Signature, LutSynthParams1::kMaxNumOutputs> func_vec;
    rinox::math::static_vector<Signature, LutSynthParams1::kMaxNumOutputs> care_vec;
    func_vec.push_back( func );
    care_vec.push_back( care );
    auto res = synth.Run(cut, func_vec, care_vec);
    if (res)
    {
      simulate<Ntk, Signature>( ntk, *res, map);
      CHECK( kitty::equal( func & care, map[res->data]&care ) );
    }
  };

  func = xs[0] & ~xs[1];
  care |= ~care;
  Check( func, care );

  func = xs[3] ^ xs[1];
  care |= ~care;
  Check( func, care );

  func = xs[0] & ~xs[1] ^ xs[3];
  care |= ~care;
  Check( func, care );

  func = xs[3] ^ xs[1] ^ xs[2];
  care |= ~care;
  Check( func, care );

  func = xs[3] ^ ((xs[1] ^ xs[2])&(xs[0]&xs[4]));
  care |= ~care;
  Check( func, care );

  func = xs[3] & ((xs[1] & xs[2])&(xs[0]&xs[4]));
  care |= ~care;
  Check( func, care );

  kitty::create_from_binary_string( func, "0001001010100001101101111011111010100011110001011110110101010011" );
  kitty::create_from_binary_string( care, "0001100011101101000011101101110110000011001000011111111110111110" );
  Check( func, care );

  for ( int cnt = 0; cnt < 100; ++cnt )
  {
    kitty::create_random( func );
    kitty::create_random( care );
    Check(func, care);
  }
}