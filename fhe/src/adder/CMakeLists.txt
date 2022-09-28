add_library(adder adder.cpp)
target_include_directories(adder PUBLIC ..)
# require at least c++17
target_compile_features(adder PUBLIC cxx_std_17)
