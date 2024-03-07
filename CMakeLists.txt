cmake_minimum_required(VERSION 3.7.1)
project(KASS
    VERSION 1.0
    DESCRIPTION "A suite of tools for k-mer analysis of RNA sequence and structure.")

# Output binaries into build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

# Build bppPipeline
add_executable(bppPipeline
    src/bppPipeline.c
    src/bppPipeline_cmdl.c
    src/rna_file_parser.c
    src/kmerHashTable.c
    src/string_utils.c
    src/utils.c
    src/parallel_helpers.c
    src/c_thread_pool/thpool.c)

# Build SKA
add_executable(SKA
    src/SKA.c
    src/rna_file_parser.c
    src/kmerHashTable.c
    src/string_utils.c
    src/utils.c
    src/parallel_helpers.c
    src/c_thread_pool/thpool.c
    src/SKA_cmdl.c)

# Include directories for SKA
target_include_directories(SKA PRIVATE src/c_thread_pool)
target_include_directories(SKA PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Include directories for BPP
target_include_directories(bppPipeline PRIVATE src/c_thread_pool)
target_include_directories(bppPipeline PRIVATE ${CMAKE_SOURCE_DIR}/include)
find_path(VIENNARNA_INCLUDE_DIRS 
    NAMES ViennaRNA/fold.h
    PATHS $ENV{HOME}/.local/include /usr/local/include)
target_include_directories(bppPipeline PRIVATE ${VIENNARNA_INCLUDE_DIRS})

# Specify compiler options based on OS
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # SKA compiler options
    target_compile_options(SKA  PRIVATE -Wall -O3 -lm)
    target_link_options(SKA  PRIVATE -Wall -O3 -flto -fopenmp -lm)

    # bppPipeline compiler options
    target_compile_options(bppPipeline PRIVATE -Wall -O3 -lm)
    target_link_options(bppPipeline PRIVATE -Wall -flto -fopenmp -O3 -lm -lstdc++)
elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin") # Apple
    target_compile_options(SKA PRIVATE -Wall -Wextra -Wpedantic -flto -O3)
    target_compile_options(bppPipeline PRIVATE -Wall -flto -O3)
else()
    message( SEND_ERROR "System ${CMAKE_SYSTEM_NAME} currently not supported.")
endif()

# Link bppPipeline with ViennaRNA library
find_library(LIBRARY_LRNA 
    NAMES RNA
    HINTS "$ENV{HOME}/.local/lib")
if(LIBRARY_LRNA)
    target_link_libraries(bppPipeline PUBLIC ${LIBRARY_LRNA})
else()
    message( SEND_ERROR "Couldn't find libRNA")
endif()

# Install files
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX /usr/local/bin)
endif()

install(TARGETS SKA bppPipeline
        DESTINATION ${CMAKE_INSTALL_PREFIX})