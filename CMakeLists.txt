cmake_minimum_required(VERSION 3.17)
project(rocketsim VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(rocketCore
    src/helper_functions.cpp
    src/physics.cpp
    src/quaternion.cpp
    src/RK4_integrator.cpp
    src/rocket_parameters.cpp
    )

target_include_directories(rocketCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/includes)

# add_executable(Simulation src/main_sim.cpp)

enable_testing()
add_executable(Test 
    tests/main_test.cpp
    tests/quaternion_test.cpp
    tests/physics_test.cpp)
target_link_libraries(Test PRIVATE rocketCore)

add_test(NAME RocketUnitTests COMMAND Test)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Environnement Debug 
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")
    # Environnement Simulation 
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
endif() 