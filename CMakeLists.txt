cmake_minimum_required(VERSION 3.7.1)
project(KATSS
    VERSION 1.0
    DESCRIPTION "A suite of tools for k-mer analysis of RNA sequence and structure.")

# Output binaries into build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

#Include zlib dependency
find_package(ZLIB REQUIRED)

# Build bppPipeline
add_executable(bppPipeline
	src/bin/bppPipeline.c
	src/bin/bppPipeline_cmdl.c
	src/rna_file_parser.c
	src/kmerHashTable.c
	src/string_utils.c
	src/utils.c
	src/parallel_helpers.c
	src/c_thread_pool/thpool.c)

# Build SKA
add_executable(SKA
	src/bin/SKA.c
	src/bin/SKA_cmdl.c
	src/rna_file_parser.c
	src/kmer_counter.c
	src/string_utils.c
	src/utils.c
	src/Regex.c)

# Build aspiire
add_executable(aspiire
	src/bin/aspiire.c
	src/bin/aspiire_cmdl.c
	src/rna_file_parser.c
	src/Regex.c
	src/ire_fold.c
	src/levenshtein.c
	src/string_utils.c
	src/utils.c
	src/parallel_helpers.c
	src/c_thread_pool/thpool.c)

# Set C standard
set_property(TARGET SKA PROPERTY C_STANDARD 11)

# Include directories for SKA
target_include_directories(SKA PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(SKA ZLIB::ZLIB)

# Include directories for BPP
target_include_directories(bppPipeline PRIVATE src/c_thread_pool)
target_include_directories(bppPipeline PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src)
find_path(VIENNARNA_INCLUDE_DIRS 
	NAMES ViennaRNA/fold.h
	PATHS $ENV{HOME}/.local/include /usr/local/include)
target_include_directories(bppPipeline PRIVATE ${VIENNARNA_INCLUDE_DIRS})

# Include directories for aspiire
target_include_directories(aspiire PRIVATE src/c_thread_pool)
target_include_directories(aspiire PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/bin)
target_include_directories(aspiire PRIVATE ${VIENNARNA_INCLUDE_DIRS})

# Specify compiler options based on OS
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	# SKA compiler options
	target_compile_options(SKA PRIVATE -Wall -O3 -lm)
	target_link_options(SKA PRIVATE -Wall -O3 -flto -fopenmp -lm)

	# bppPipeline compiler options
	target_compile_options(bppPipeline PRIVATE -Wall -O3 -lm)
	target_link_options(bppPipeline PRIVATE -Wall -flto -fopenmp -O3 -lm -lstdc++)

	# aspiire compiler options
	target_compile_options(aspiire PRIVATE -Wall -Wextra -O3)
	target_link_options(aspiire PRIVATE -Wall -Wextra -flto -fopenmp -O3 -lm -lstdc++)
elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin") # Apple
	target_compile_options(SKA PRIVATE -Wall -Wextra -Wpedantic -flto -O3)
	target_compile_options(bppPipeline PRIVATE -Wall -flto -O3)
	target_compile_options(aspiire PRIVATE -Wall -Wextra -flto -O3)
else()
	message( SEND_ERROR "System ${CMAKE_SYSTEM_NAME} currently not supported.")
endif()

# Link bppPipeline & aspiire with libraries
find_library(LIBRARY_LRNA 
	NAMES RNA
	HINTS "$ENV{HOME}/.local/lib")
if(LIBRARY_LRNA)
	target_link_libraries(bppPipeline PUBLIC ${LIBRARY_LRNA} ZLIB::ZLIB)
	target_link_libraries(aspiire PUBLIC ${LIBRARY_LRNA} ZLIB::ZLIB)
else()
	message( SEND_ERROR "Couldn't find libRNA")
endif()

# Install files
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX /usr/local/bin)
endif()

install(TARGETS SKA bppPipeline aspiire
        DESTINATION ${CMAKE_INSTALL_PREFIX})
