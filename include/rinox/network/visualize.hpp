#pragma once
#include <drag/graph.hpp>
#include <drag/layout.hpp>
#include <fstream>
#include <limits>
#include <mockturtle/mockturtle.hpp> // for foreach_fanin etc.
#include <string>
#include <unordered_map>

namespace rinox::graphics
{

struct NodeStyle
{
  std::string fill = "lightblue";
  std::string stroke = "black";
};

struct NodeOutput
{
  uint64_t node;
  uint64_t output_index;
  bool operator==( const NodeOutput& o ) const noexcept
  {
    return node == o.node && output_index == o.output_index;
  }
};

struct NodeOutputHash
{
  std::size_t operator()( const NodeOutput& no ) const noexcept
  {
    return ( static_cast<std::size_t>( no.node ) << 8 ) ^ no.output_index;
  }
};

template<typename Ntk, typename Window>
inline void visualize_multioutput_window( Ntk const& ntk, Window const& w, std::string const& filename )
{
  drag::graph g;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  std::unordered_map<NodeOutput, uint32_t, NodeOutputHash> id_map;
  std::unordered_map<uint32_t, std::string> labels;
  std::unordered_map<uint32_t, NodeStyle> styles;

  auto register_signal_outputs = [&]( signal const& f, std::string const& color ) {
    node n = ntk.get_node( f );
    NodeOutput key{ n, f.output };
    if ( !id_map.count( key ) )
    {
      uint32_t node_id = g.add_node();
      id_map[key] = node_id;
      labels[node_id] = std::to_string( f.index ) + "_" + std::to_string( f.output );
      styles[node_id].fill = color;
    }
  };

  auto register_node_outputs = [&]( node const& n, std::string const& color ) {
    ntk.foreach_output( n, [&]( auto const& f ) {
      NodeOutput key{ n, f.output };
      if ( !id_map.count( key ) )
      {
        uint32_t node_id = g.add_node();
        id_map[key] = node_id;
        labels[node_id] = std::to_string( f.index ) + "_" + std::to_string( f.output );
        styles[node_id].fill = color;
      }
    } );
  };

  // Assign colors per category
  for ( auto const& f : w.inputs )
    register_signal_outputs( f, "#a1d99b" ); // green inputs
  for ( auto const& f : w.divs )
    register_signal_outputs( f, "#9ecae1" );   // blue divisors
  register_node_outputs( w.pivot, "#e6550d" ); // pivot = red/orange
  for ( auto const& n : w.mffc )
    register_node_outputs( n, "#fc9272" ); // red mffc
  for ( auto const& n : w.tfos )
    register_node_outputs( n, "#fdae6b" ); // orange tfo
  for ( auto const& f : w.outputs )
    register_signal_outputs( f, "#756bb1" ); // purple outputs

  // Take a snapshot of the registered outputs
  std::vector<NodeOutput> keys;
  keys.reserve( id_map.size() );
  for ( auto const& [key, _] : id_map )
  {
    keys.push_back( key );
  }

  for ( auto const& key : keys )
  {
    auto dst = id_map.at( key );
    auto const n = key.node;

    ntk.foreach_output( n, [&]( auto const& f ) {
      if ( f.output != key.output_index )
        return;

      // Fanins
      ntk.foreach_fanin( f, [&]( auto const& fi, auto ) {
        NodeOutput src_key{ ntk.get_node( fi ), fi.output };
        // if (!id_map.count(src_key)) {
        //   uint32_t node_id = g.add_node();
        //   id_map[src_key] = node_id;
        //   labels[node_id] = std::to_string(fi.index) + "_" + std::to_string(fi.output);
        //   styles[node_id].fill = "#ff8fe9";
        // }
        // g.add_edge(id_map[src_key], dst);

        if ( id_map.count( src_key ) )
          g.add_edge( id_map[src_key], dst );
      } );

      // Fanouts
      ntk.foreach_fanout( f, [&]( auto const& no ) {
        ntk.foreach_output( no, [&]( auto const& fo ) {
          NodeOutput fanout_key{ no, fo.output };
          // if (!id_map.count(fanout_key)) {
          //   uint32_t node_id = g.add_node();
          //   id_map[fanout_key] = node_id;
          //   labels[node_id] = std::to_string(fo.index) + "_" + std::to_string(fo.output);
          //   styles[node_id].fill = "#ff8fe9";
          // }
          // g.add_edge(dst, id_map[fanout_key]);
          if ( id_map.count( fanout_key ) )
            g.add_edge( dst, id_map[fanout_key] );
        } );
      } );
    } );
  }

  // Compute Sugiyama layout
  drag::sugiyama_layout layout( g );
  auto const& verts = layout.vertices();

  // Compute bounding box
  double min_x = std::numeric_limits<double>::max();
  double min_y = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double max_y = std::numeric_limits<double>::lowest();

  for ( auto const& v : verts )
  {
    min_x = std::min( min_x, static_cast<double>( v.pos.x ) );
    min_y = std::min( min_y, static_cast<double>( v.pos.y ) );
    max_x = std::max( max_x, static_cast<double>( v.pos.x ) );
    max_y = std::max( max_y, static_cast<double>( v.pos.y ) );
  }

  double margin = 50.0;
  min_x -= margin;
  min_y -= margin;
  max_x += margin;
  max_y += margin;

  // Write SVG
  std::ofstream os( filename + ".svg" );
  os << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" "
     << "viewBox=\"" << min_x << " " << min_y << " "
     << ( max_x - min_x ) << " " << ( max_y - min_y ) << "\" "
     << "width=\"" << ( max_x - min_x ) << "\" "
     << "height=\"" << ( max_y - min_y ) << "\">\n";

  // Edges
  for ( auto const& path : layout.edges() )
  {
    auto const& pts = path.points;
    for ( size_t i = 1; i < pts.size(); ++i )
    {
      os << "<line x1=\"" << pts[i - 1].x << "\" y1=\"" << pts[i - 1].y
         << "\" x2=\"" << pts[i].x << "\" y2=\"" << pts[i].y
         << "\" stroke=\"black\" stroke-width=\"1\" />\n";
    }
  }

  // Nodes
  for ( auto const& v : verts )
  {
    uint32_t id = v.u;
    auto lit = labels.find( id );
    auto sit = styles.find( id );
    std::string full_label = ( lit != labels.end() ) ? lit->second : std::to_string( id );
    std::string short_label = ( full_label.size() > 6 ) ? full_label.substr( 0, 6 ) + "…" : full_label;
    NodeStyle style = ( sit != styles.end() ) ? sit->second : NodeStyle{};

    os << "<circle cx=\"" << v.pos.x << "\" cy=\"" << v.pos.y
       << "\" r=\"12\" fill=\"" << style.fill
       << "\" stroke=\"" << style.stroke << "\">\n";
    os << "  <title>" << full_label << "</title>\n";
    os << "</circle>\n";

    os << "<text x=\"" << v.pos.x << "\" y=\"" << v.pos.y
       << "\" font-size=\"9\" text-anchor=\"middle\" dy=\".35em\">"
       << short_label << "</text>\n";
  }

  os << "</svg>\n";
}

} // namespace rinox::graphics
