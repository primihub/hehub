add_library(common rnspolynomial.cpp bigint.cpp bigintpoly.cpp mod_arith.cpp ntt.cpp)
# require at least c++17
target_compile_features(common PUBLIC cxx_std_17)