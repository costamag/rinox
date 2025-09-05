#pragma once
#include <drag/graph.hpp>
#include <drag/layout.hpp>
#include <drag/svg_writer.hpp>
#include <string>

namespace rinox::graphics {

/**
 * @brief Visualize a DAG into an SVG file using Drag's automatic layout
 *
 * @param g        The DAG to render
 * @param filename Path to the output SVG file
 */
inline void visualize_dag(const drag::Graph& g, const std::string& filename)
{
    drag::Layout layout(g);
    layout.run(); // automatic Sugiyama DAG layout
    drag::SvgWriter::write(filename, g, layout.positions());
}

} // namespace rinox::graphics
