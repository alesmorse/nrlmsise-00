cmake_minimum_required (VERSION 3.7.0 FATAL_ERROR)

project(nrlmsise00
    LANGUAGES CXX
    VERSION 1.0.0)

# Compiler options
set(CMAKE_CXX_FLAGS "-Wall -Wextra")

message(STATUS "Project will be installed to ${CMAKE_INSTALL_PREFIX}")

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type is set to ${CMAKE_BUILD_TYPE}")

# General testing information
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
include(CTest)

add_library(${PROJECT_NAME}
    SHARED
        include/Nrlmsise00.hpp
        src/Nrlmsise00.cpp
    )

set_target_properties(${PROJECT_NAME}
    PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
        POSITION_INDEPENDENT_CODE 1
    )

target_sources(${PROJECT_NAME}
    PRIVATE
        src/Nrlmsise00_p.cpp
        src/Nrlmsise00_math.cpp
        src/Nrlmsise00_data.cpp
    )

target_include_directories(${PROJECT_NAME}
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/detail>
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
    )

install(TARGETS  ${PROJECT_NAME}
        LIBRARY  DESTINATION lib
        ARCHIVE  DESTINATION lib
        INCLUDES DESTINATION include/${PROJECT_NAME})

install(DIRECTORY include/ DESTINATION include/${PROJECT_NAME}
    FILES_MATCHING PATTERN "*.hpp")

# Unit tests
enable_testing()
add_subdirectory(tests)

# Doxygen documentation
find_package(Doxygen)
if (DOXYGEN_FOUND)
    include(UseDoxygenDoc)
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
endif()