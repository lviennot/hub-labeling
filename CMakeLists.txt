cmake_minimum_required(VERSION 3.0)

project(contraction_hierarchies CXX)

include(CMakeDetermineCXXCompiler)

# default to RelWithDebInfo
# to debug : cmake -DCMAKE_BUILD_TYPE:STRING=Debug ..
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
	set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

message(STATUS "  Build type ${CMAKE_BUILD_TYPE}")
message(STATUS "  Compiler ID : ${CMAKE_CXX_COMPILER_ID}")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(EXTRA_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -pedantic" CACHE STRING "Extra flags used by the compiler during all build types.")
        message(STATUS "  Extra flags : ${EXTRA_CXX_FLAGS}")
	set(EXTRA_EXE_LINKER_FLAGS "" CACHE STRING "Extra flags used by the linker.")
	set(EXTRA_EXE_LINKER_FLAGS_RELEASE "-flto" CACHE STRING "Extra flags used by the linker for the Release build type.")
	set(EXTRA_EXE_LINKER_FLAGS_RELWITHDEBINFO "-flto" CACHE STRING "Extra flags used by the linker for the RelWithDebInfo build type.")
endif()

# 
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -DNDEBUG")
message(STATUS "  Flags RelWithDebInfo: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS} -std=c++17 -fpermissive")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${EXTRA_EXE_LINKER_FLAGS} -pthread")

message(STATUS "  CXX Flags: ${CMAKE_CXX_FLAGS}")


add_executable(hl_trans
        src/hl_trans.cc
)

add_executable(unit
        src/heap_unit.cc
        src/mgraph_unit.cc
        src/unit.cc
)




