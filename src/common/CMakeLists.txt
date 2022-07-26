add_library(common rnspolynomial.cpp bigint.cpp bigintpoly.cpp)
# require at least c++17
target_compile_features(adder PUBLIC cxx_std_17)