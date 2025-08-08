cmake_minimum_required(VERSION 3.10.0)
project(testC VERSION 0.1.0 LANGUAGES C CXX)

add_executable(stateMachine)
target_include_directories(stateMachine PRIVATE inc)
target_sources(stateMachine PRIVATE src/stateMachine.cpp src/main.cpp)

add_executable(bisection)
target_include_directories(bisection PRIVATE inc)
target_sources(
    bisection 
    PRIVATE 
    src/bisection.c 
    src/drvBisection.c 
    src/bisectionMain.c
    src/io.c
    src/error.c
    src/operator.c
    src/convert.c
)

add_executable(qzexplicit)
target_include_directories(qzexplicit PRIVATE inc)
target_sources(
  qzexplicit
  PRIVATE
  src/drvEigen.c
  src/qzExplicit.c
  src/qzExplicitMain.c
  src/io.c
  src/error.c
  src/operator.c
  src/convert.c
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  # using Clang
  target_link_libraries(bisection m)
  target_link_libraries(qzexplicit m)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # using GCC
  target_link_libraries(bisection m)
  target_link_libraries(qzexplicit m)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  # using Intel C++
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # using Visual Studio C++
endif()
message(">>>> Using Compiler: ${CMAKE_CXX_COMPILER_ID} <<<<")

