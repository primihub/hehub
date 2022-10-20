add_library(${PROJECT_NAME})
# require at least c++17
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

set(THIRD_PARTY_DIR ${PROJECT_SOURCE_DIR}/third-party)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${THIRD_PARTY_DIR}/range-v3)

Option(HEHUB_DEBUG_FHE OFF)
if(HEHUB_DEBUG)
    set(HEHUB_DEBUG_FHE ON)
endif()
if(HEHUB_DEBUG_FHE)
    target_compile_definitions(${PROJECT_NAME} PUBLIC HEHUB_DEBUG_FHE)
endif()

Option(HEHUB_DEBUG_RLWE_ZERO_C1 OFF)
Option(HEHUB_DEBUG_RLWE_ZERO_E OFF)
if(HEHUB_DEBUG_RLWE_ZERO_C1)
    target_compile_definitions(${PROJECT_NAME} PUBLIC HEHUB_DEBUG_RLWE_ZERO_C1)
endif()
if(HEHUB_DEBUG_RLWE_ZERO_E)
    target_compile_definitions(${PROJECT_NAME} PUBLIC HEHUB_DEBUG_RLWE_ZERO_E)
endif()

add_subdirectory(common)
add_subdirectory(primitives)
add_subdirectory(bgv)
add_subdirectory(ckks)
