find_package(pybind11 REQUIRED)

pybind11_add_module(nrlmsise00_py bindings.cpp)

target_link_libraries(nrlmsise00_py PRIVATE nrlmsise00)

set_target_properties(nrlmsise00_py
    PROPERTIES
        OUTPUT_NAME nrlmsise00
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
    )

install(TARGETS nrlmsise00_py DESTINATION .)
