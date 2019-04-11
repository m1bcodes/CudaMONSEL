project(CudaMONSEL
VERSION 1.0
DESCRIPTION "CUDA MONSLE"
LANGUAGES CXX CUDA
)

message(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(CudaMONSEL/gov/nist/microanalysis)
