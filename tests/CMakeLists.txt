FIND_PACKAGE(Cppunit REQUIRED)

list(APPEND test_FILES
    Nrlmsise00Test.cpp)

add_executable(unit_tests
    TestRunner.cpp
    ${test_FILES}
    )

TARGET_LINK_LIBRARIES(unit_tests
    ${CPPUNIT_LIBRARIES}
    nrlmsise00
    )

add_test(Nrlmsise00 unit_tests)