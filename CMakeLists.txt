cmake_minimum_required(VERSION 3.10.0)
project(testC VERSION 0.1.0 LANGUAGES C CXX)

add_executable(stateMachine)
target_include_directories(stateMachine PRIVATE inc)
target_sources(stateMachine PRIVATE src/stateMachine.cpp src/main.cpp)

