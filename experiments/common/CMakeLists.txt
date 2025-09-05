file(GLOB_RECURSE RINOX_COMMON_SOURCES CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_library(rinox_exp_common STATIC ${RINOX_COMMON_SOURCES})

target_include_directories(rinox_exp_common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/lib/mockturtle  # <-- ensure mockturtle/ is found
)

target_link_libraries(rinox_exp_common PUBLIC
  rinox_inspect
  mockturtle
  rinox_headers
  rinox_build_options
)