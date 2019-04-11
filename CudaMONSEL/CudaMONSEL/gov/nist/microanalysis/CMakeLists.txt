project(Microanalysis
VERSION 1.0 
DESCRIPTION "Microanalysis"
LANGUAGES CXX CUDA
)

include(CTest)

set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(UTILITY ${CMAKE_CURRENT_SOURCE_DIR}/Utility)
set(EPQTESTS ${CMAKE_CURRENT_SOURCE_DIR}/EPQTests)
set(AMPHIBIAN_ROOT ${CMAKE_SOURCE_DIR}/Amphibian/Amphibian/Amphibian)

message(${PROJECT_ROOT})
message(${EPQTESTS})
message(${UTILITY})

string(APPEND NVCC_FLAGS " -keep")

# GLOB
file(GLOB UTILITY_UNCERTAINVALUE2 ${UTILITY}/UncertainValue2.cu)

# lib
add_library(utility-uncertainvalue2 ${UTILITY_UNCERTAINVALUE2})
target_include_directories(utility-uncertainvalue2 PUBLIC ${AMPHIBIAN_ROOT})
target_include_directories(utility-uncertainvalue2 PUBLIC ${UTILITY})
set_target_properties(utility-uncertainvalue2 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_compile_features(utility-uncertainvalue2 PUBLIC cxx_std_11)

# tests
if(BUILD_TESTING)
    # file(GLOB UTILITY_TEST_UNCERTAINVALUE2 ${EPQTESTS}/UncertainValue2Test.cu)

    # add_library(utility-test-uncertainvalue2 ${UTILITY_TEST_UNCERTAINVALUE2})
    # target_include_directories(utility-test-uncertainvalue2 PUBLIC ${AMPHIBIAN_ROOT})
    # target_include_directories(utility-test-uncertainvalue2 PUBLIC ${UTILITY})
    # target_include_directories(utility-test-uncertainvalue2 PUBLIC ${EPQTESTS})
    # set_target_properties(utility-test-uncertainvalue2 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    # target_compile_features(utility-test-uncertainvalue2 PUBLIC cxx_std_11)

    # exe
    # add_executable(amphibian ${MAIN_SRCS})
    # set_target_properties(amphibian PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #target_link_libraries(amphibian amphibian-lib)
endif()
