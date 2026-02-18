find_package(pybind11 REQUIRED)

pybind11_add_module(nrlmsise00_py
    bindings.cpp
    ${CMAKE_SOURCE_DIR}/src/Nrlmsise00.cpp
    ${CMAKE_SOURCE_DIR}/src/Nrlmsise00_p.cpp
    ${CMAKE_SOURCE_DIR}/src/Nrlmsise00_math.cpp
    ${CMAKE_SOURCE_DIR}/src/Nrlmsise00_data.cpp
)

target_include_directories(nrlmsise00_py PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/detail
)

set_target_properties(nrlmsise00_py
    PROPERTIES
        OUTPUT_NAME nrlmsise00
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
    )

install(TARGETS nrlmsise00_py DESTINATION .)
