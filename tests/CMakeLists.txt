add_executable(tests tests.cpp adder_t.cpp common.cpp bigint_t.cpp mod_arith_t.cpp)
target_link_libraries(tests PUBLIC adder)
target_link_libraries(tests PUBLIC common)
target_include_directories(tests PUBLIC ${PROJECT_SOURCE_DIR}/third-party)

# allow user to run tests with `make test` or `ctest`
include(../cmake/Catch.cmake)
catch_discover_tests(tests)
