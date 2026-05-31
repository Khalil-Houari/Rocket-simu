cmake_minimum_required(VERSION 3.10)
project(rocketsim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(rocketsim main.cpp)

# Check which compiler is being used and apply the right flags
if(MSVC)
    # /W4 is the MSVC equivalent to GCC's -Wall -Wextra
    target_compile_options(rocketsim PRIVATE /W4)
else()
    # Link math library on Unix/Linux systems (not needed on MSVC)
    target_link_libraries(rocketsim m)
endif()