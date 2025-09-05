# Add common lib
add_subdirectory(common)

# Build each experiment .cpp (excluding common/)
file(GLOB_RECURSE EXPERIMENT_CPP CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

foreach(SRC ${EXPERIMENT_CPP})
  # Skip common/
  string(FIND "${SRC}" "/common/" COMMON_POS)
  if (NOT COMMON_POS EQUAL -1)
    continue()
  endif()

  get_filename_component(NAME_WE ${SRC} NAME_WE)
  file(RELATIVE_PATH REL_DIR ${CMAKE_CURRENT_SOURCE_DIR} ${SRC})
  get_filename_component(DIR ${REL_DIR} DIRECTORY)
  set(OUTDIR ${CMAKE_BINARY_DIR}/experiments/${DIR})

  add_executable(${NAME_WE} ${SRC})
  set_target_properties(${NAME_WE} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${OUTDIR}
  )

  target_include_directories(${NAME_WE} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_SOURCE_DIR}/src
  )

target_link_libraries(rinox_exp_common PUBLIC mockturtle)

  target_link_libraries(${NAME_WE} PRIVATE
    rinox_inspect
    rho_cli
    rinox_exp_common
    rinox_build_options
  )
endforeach()
