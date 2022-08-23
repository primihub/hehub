project(hehub)
add_library(${PROJECT_NAME})
# require at least c++17
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(adder)
add_subdirectory(common)
add_subdirectory(primitives)
