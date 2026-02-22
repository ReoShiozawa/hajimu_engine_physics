cmake_minimum_required(VERSION 3.20)
project(jp_engine_physics C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ─── hajimu_plugin.h ───────────────────────────────────
find_path(HAJIMU_INCLUDE_DIR hajimu_plugin.h
    PATHS
        ${CMAKE_SOURCE_DIR}/../../jp/include
        ${CMAKE_SOURCE_DIR}/../jp/include
        /usr/local/include/hajimu
    REQUIRED)
message(STATUS "hajimu include: ${HAJIMU_INCLUDE_DIR}")

# ─── ターゲット ────────────────────────────────────────
set(SRCS
    src/phys.c
    src/plugin.c
)

add_library(engine_physics SHARED ${SRCS})

set_target_properties(engine_physics PROPERTIES
    PREFIX  ""
    SUFFIX  ".hjp"
    OUTPUT_NAME "engine_physics"
)

target_include_directories(engine_physics PRIVATE
    ${HAJIMU_INCLUDE_DIR}
    include
    src
)

# ─── リンク ────────────────────────────────────────────
if (APPLE)
    target_compile_definitions(engine_physics PRIVATE _USE_MATH_DEFINES)
    target_link_libraries(engine_physics PRIVATE
        "-framework CoreFoundation"
        m
    )
else()
    target_link_libraries(engine_physics PRIVATE m)
endif()

target_compile_options(engine_physics PRIVATE
    -O2 -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers
)
