cmake_minimum_required(VERSION 3.14)

# Set a name and a version number for your project:
project(hehub VERSION 0.1 LANGUAGES CXX)

# this needs to be in the top level CMakeLists.txt to enable tests
include(CTest)

set(CMAKE_VERBOSE_MAKEFILE OFF)

# ── Architecture-native optimizations ──────────────────────────────────
option(HEHUB_NATIVE_ARCH "Enable -march=native for CPU-specific optimizations" ON)
if(HEHUB_NATIVE_ARCH)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-march=native" HAS_MARCH_NATIVE)
    if(HAS_MARCH_NATIVE)
        add_compile_options(-march=native)
    endif()
endif()

# ── Interprocedural optimization / LTO ─────────────────────────────────
option(HEHUB_LTO "Enable link-time optimization" ON)
if(HEHUB_LTO)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-flto" HAS_LTO)
    if(HAS_LTO)
        add_compile_options(-flto)
        add_link_options(-flto)
    endif()
endif()

# ── OpenMP ─────────────────────────────────────────────────────────────
if(APPLE)
    # Homebrew libomp (keg-only) on Apple Silicon
    set(OPENMP_PREFIX "/opt/homebrew/opt/libomp")
    if(EXISTS "${OPENMP_PREFIX}")
        set(OpenMP_CXX_FLAGS "-Xclang -fopenmp -I${OPENMP_PREFIX}/include")
        set(OpenMP_CXX_LIB_NAMES "libomp")
        set(OpenMP_libomp_LIBRARY "${OPENMP_PREFIX}/lib/libomp.dylib")
    endif()
endif()
find_package(OpenMP QUIET)

# ── Debug mode ─────────────────────────────────────────────────────────
Option(HEHUB_DEBUG OFF)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
    set(HEHUB_DEBUG ON)
endif()
if(HEHUB_DEBUG)
    add_compile_definitions(HEHUB_DEBUG)
endif()

# compile the libraries
add_subdirectory(src)

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
