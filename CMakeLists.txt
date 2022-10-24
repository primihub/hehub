cmake_minimum_required(VERSION 3.14)

# Set a name and a version number for your project:
project(hehub VERSION 0.1 LANGUAGES CXX)

# this needs to be in the top level CMakeLists.txt to enable tests
include(CTest)

set(CMAKE_VERBOSE_MAKEFILE OFF)

Option(HEHUB_DEBUG OFF)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
    set(HEHUB_DEBUG ON)
endif()
if(HEHUB_DEBUG)
    add_compile_definitions(HEHUB_DEBUG)
endif()

# compile the FHE library
add_subdirectory(src/fhe)

# compile the tests
add_subdirectory(tests)

# compile the benchmarks
add_subdirectory(bench)

# compile the examples
add_subdirectory(examples)

# optionally add doxygen target to generate documentation
option(BUILD_DOCS "Enable building of documentation (requires Doxygen)" OFF)
if(BUILD_DOCS)
	find_package(Doxygen REQUIRED)
	set(DOXYGEN_EXCLUDE_PATTERNS "${CMAKE_SOURCE_DIR}/third-party/*")
	doxygen_add_docs(doxygen
		${CMAKE_SOURCE_DIR}
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
endif()