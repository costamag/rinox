cmake_minimum_required(VERSION 3.16)

# Collect any potential .cpp sources in diagnostics + visualization
file(GLOB_RECURSE RINOX_INSPECT_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/rinox/diagnostics/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/rinox/graphics/*.cpp"
)

set(RINOX_INSPECT_PUBLIC_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# If we have no sources, make it header-only
if (RINOX_INSPECT_SOURCES STREQUAL "")
    add_library(rinox_inspect INTERFACE)
    target_include_directories(rinox_inspect INTERFACE
        $<BUILD_INTERFACE:${RINOX_INSPECT_PUBLIC_DIR}>
        $<INSTALL_INTERFACE:include>
    )
    target_compile_features(rinox_inspect INTERFACE cxx_std_17)
else()
    add_library(rinox_inspect STATIC ${RINOX_INSPECT_SOURCES})
    target_include_directories(rinox_inspect PUBLIC
        $<BUILD_INTERFACE:${RINOX_INSPECT_PUBLIC_DIR}>
        $<INSTALL_INTERFACE:include>
    )
    target_compile_features(rinox_inspect PUBLIC cxx_std_17)
endif()

# Link dependencies for diagnostics/visualization features
if (TARGET mockturtle::mockturtle)
    target_link_libraries(rinox_inspect INTERFACE mockturtle::mockturtle)
elseif (TARGET mockturtle)
    target_link_libraries(rinox_inspect INTERFACE mockturtle)
endif()

if (TARGET lorina::lorina)
    target_link_libraries(rinox_inspect INTERFACE lorina::lorina)
endif()

if (TARGET fmt::fmt)
    target_link_libraries(rinox_inspect INTERFACE fmt::fmt)
endif()