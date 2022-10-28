add_library(${PROJECT_NAME}-circuits SHARED linear_algebra.cpp)
# require at least c++17
target_compile_features(${PROJECT_NAME}-circuits PUBLIC cxx_std_17)

set(THIRD_PARTY_DIR ${PROJECT_SOURCE_DIR}/third-party)

target_include_directories(${PROJECT_NAME}-circuits PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME}-circuits PUBLIC ${THIRD_PARTY_DIR}/range-v3)
target_link_libraries(${PROJECT_NAME}-circuits PUBLIC ${PROJECT_NAME})

