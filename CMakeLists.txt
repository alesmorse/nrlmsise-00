cmake_minimum_required (VERSION 3.7.0 FATAL_ERROR)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)


project(NRLMSISE-00 LANGUAGES CXX)
set(LIB_NAME "nrlmsise00")

# General testing information
include(CTest)
include(UseDoxygenDoc)

# Compiler options
set(CMAKE_CXX_FLAGS "-Wall -Wextra")

add_library(${LIB_NAME}
    SHARED
        include/Nrlmsise00.hpp
        src/Nrlmsise00.cpp
    )

set_target_properties(${LIB_NAME}
    PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
        POSITION_INDEPENDENT_CODE 1
    )

target_sources(${LIB_NAME}
    PRIVATE
        src/Nrlmsise00_p.cpp
        src/Nrlmsise00_math.cpp
        src/Nrlmsise00_data.cpp
    )

target_include_directories(${LIB_NAME}
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/detail>
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/${LIB_NAME}>
    )

install(TARGETS  ${LIB_NAME}
        EXPORT   atmosTargets
        LIBRARY  DESTINATION lib
        ARCHIVE  DESTINATION lib
        INCLUDES DESTINATION include/${LIB_NAME})

install(DIRECTORY include/ DESTINATION include/${LIB_NAME}
    FILES_MATCHING PATTERN "*.hpp")

add_subdirectory(tests)

add_doxygen_doc(
    BUILD_DIR
        ${CMAKE_CURRENT_BINARY_DIR}/docs
    DOXY_FILE
        ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in
    TARGET_NAME
        docs
    COMMENT
        "HTML Documentation"
)