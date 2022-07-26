cmake_minimum_required(VERSION 3.8.2)

# Set a name and a version number for your project:
project(cpp-project-template VERSION 0.0.1 LANGUAGES CXX)

# this needs to be in the top level CMakeLists.txt to enable tests
include(CTest)

# compile the library
add_subdirectory(src)

# compile the tests
add_subdirectory(tests)

# compile the application
add_subdirectory(app)

# optionally add doxygen target to generate documentation
option(BUILD_DOCS "Enable building of documentation (requires Doxygen)" OFF)
if(BUILD_DOCS)
	find_package(Doxygen REQUIRED)
	set(DOXYGEN_EXCLUDE_PATTERNS "${CMAKE_SOURCE_DIR}/ext/*")
	doxygen_add_docs(doxygen
		${CMAKE_SOURCE_DIR}
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
endif()