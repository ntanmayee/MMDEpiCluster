cmake_minimum_required(VERSION 3.12)
project(mmdiff3)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE 1)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

message(STATUS "The build directory is ${CMAKE_SOURCE_DIR}")

find_package(OPENMP)
if(OpenMP_FOUND)
    message(STATUS "OpenMP Found")

    set(OpenMP_CXX_FLAGS "-Xclang -fopenmp")
    set(OpenMP_C_FLAGS "-Xclang -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

    if(NOT TARGET OpenMP::OpenMP_CXX)
        add_library(OpenMP_TARGET INTERFACE)
        add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
        target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
        find_package(Threads REQUIRED)
        target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
        target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
    endif()
    list(APPEND LINK_LIBRARIES OpenMP::OpenMP_CXX)
else()
    message(STATUS "Failed to locate OpenMP, proceeding without")
endif()

set(SOURCES
        mmd.hpp
        mmd.cpp
        kernel_function.hpp
        rbf_joint_kernel.hpp
        rbf_joint_kernel.cpp
        rbf_joint_discrete_kernel.hpp
        rbf_joint_discrete_kernel.cpp
        mmdiff3.hpp
        mmdiff3.cpp
        R_init.hpp
        R_init.cpp)

find_package(LibR REQUIRED)
message(STATUS "Including Rlib directories ${LIBR_INCLUDE_DIRS}")
include_directories( PRIVATE ${LIBR_INCLUDE_DIRS})

add_library(mmdiff3 SHARED ${SOURCES})
target_link_libraries(mmdiff3 PRIVATE ${LIBR_LIBRARIES} ${LINK_LIBRARIES})

set_target_properties(mmdiff3 PROPERTIES COMPILE_FLAGS "-ggdb")
set_target_properties(mmdiff3 PROPERTIES PREFIX "")

if(NOT WIN32) # Need to force .so on Mac, but Windows needs to be .dll
    set_target_properties(mmdiff3 PROPERTIES SUFFIX ".so")
    install(TARGETS mmdiff3 LIBRARY DESTINATION . )
else(NOT WIN32)
    install(TARGETS mmdiff3 RUNTIME DESTINATION . )
endif(NOT WIN32)