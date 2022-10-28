add_executable(tests tests.cpp common_t.cpp bigint_t.cpp 
    mod_arith_t.cpp ntt_t.cpp rlwe_t.cpp bgv_t.cpp ckks_t.cpp lin_alg_t.cpp)
target_link_libraries(tests PUBLIC hehub)
target_link_libraries(tests PUBLIC hehub-circuits)
target_include_directories(tests PUBLIC ${PROJECT_SOURCE_DIR}/third-party)

# allow user to run tests with `make test` or `ctest`
include(../cmake/Catch.cmake)
catch_discover_tests(tests)
