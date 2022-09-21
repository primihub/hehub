target_sources(${PROJECT_NAME} PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/rlwe.cpp)

add_subdirectory(bgv)
add_subdirectory(ckks)
