project(SALT_PROJECT)
cmake_minimum_required(VERSION 3.5.1)

# CMake Options
set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_COMPILER ${CMAKE_CXX_COMPILER}) # use g++ for .c files (flute)


set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
set(CMAKE_BUILD_TYPE "Release")
# set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")
# set(CMAKE_BUILD_TYPE "Debug")


set(POWV9_DAT ${PROJECT_SOURCE_DIR}/src/salt/base/flute/POWV9.dat)
set(POST9_DAT ${PROJECT_SOURCE_DIR}/src/salt/base/flute/POST9.dat)

include_directories(${PROJECT_SOURCE_DIR}/src)
file(GLOB SALT_SRCS src/salt/*.cpp 
                    src/salt/*/*.cpp 
                    src/salt/base/flute/*.c 
                    src/salt/refine/*.cpp 
                    src/other_methods/*.cpp
                    src/delay_tree/*.cpp)
file(GLOB OTHER_METHODS_SRCS src/other_methods/*.cpp)
add_library(salt STATIC ${SALT_SRCS})
add_library(other_methods OBJECT ${OTHER_METHODS_SRCS})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../run)

# eval_single_salt
add_executable(eval_single_salt src/eval/eval_single_main.cpp $<TARGET_OBJECTS:other_methods>)
target_link_libraries(eval_single_salt salt)

# eval_batch_salt
add_executable(eval_batch_salt src/eval/eval_batch_main.cpp $<TARGET_OBJECTS:other_methods>)
target_link_libraries(eval_batch_salt salt)

#run
add_executable(eval_specified_nets src/eval/eval_specified_nets.cpp $<TARGET_OBJECTS:other_methods>)
target_link_libraries(eval_specified_nets salt)

#benchmark
add_executable(benchmark_constuct src/benchmark/benchmark_constuct.cpp $<TARGET_OBJECTS:other_methods>)
target_link_libraries(benchmark_constuct salt)