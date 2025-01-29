cmake_minimum_required(VERSION 3.5.0)
project(BFES VERSION 0.1.0 LANGUAGES C CXX)

# Unless specified otherwise, build for release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Link time optimization
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-flto" COMPILER_SUPPORTS_LTO)
if(COMPILER_SUPPORTS_LTO)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -flto")
else()
    message(WARNING "LTO is not supported by the compiler")
endif()

# Macro for fetching all header files
MACRO(HEADER_DIRECTORIES return_list)
    FILE(GLOB_RECURSE new_list "*.hpp" "*.h")
    SET(dir_list "")
    FOREACH(file_path ${new_list})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()

# Require C++ 20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Basic flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -fno-finite-math-only")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -march=native -fno-finite-math-only")

# Special OMP Flags
if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xclang -fopenmp")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
endif()

# FIND OPENMP
if(APPLE)
    set(LIBOMP_DIR "/opt/homebrew/Cellar/libomp/19.1.5")
    if (NOT LIBOMP_INCLUDE_DIR)
        set(LIBOMP_INCLUDE_DIR "${LIBOMP_DIR}/include")
    endif()
    if (NOT LIBOMP_LIB)
        set(LIBOMP_LIB "${LIBOMP_DIR}/lib/libomp.dylib") 
    endif()
    find_path(LIBOMP_INCLUDE_DIR omp.h REQUIRED)
    find_library(LIBOMP_LIB libomp REQUIRED)
    if (LIBOMP_INCLUDE_DIR AND LIBOMP_LIB)
        message(STATUS "OpenMP found")
        set(LIBOMP_LIBS ${LIBOMP_LIB})
    else()
        message(FATAL_ERROR "OpenMP library not found")
    endif()
else()
    find_package(OpenMP REQUIRED)
    if (OpenMP_CXX_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        message(STATUS "OpenMP found")
        set(LIBOMP_INCLUDE_DIR "")
    else()
        message(FATAL_ERROR "OpenMP library not found")
    endif()
endif()

# FIND MOSEK
set(MOSEK_DIR "/home/francesco-dadalt/Documents/ETH/Coding/MOSEK/mosek/10.2/tools/platform/linux64x86")
# set(MOSEK_DIR "/opt/mosek/10.0/tools/platform/osxaarch64")
set(MOSEK_INCLUDE_DIR "${MOSEK_DIR}/h/")
if (APPLE)
    set(MOSEK_LIB "${MOSEK_DIR}/bin/libfusion64.dylib")
else()
    set(MOSEK_LIB "${MOSEK_DIR}/bin/libfusion64.so")
endif()
find_path(MOSEK_INCLUDE_DIR mosek.h REQUIRED)
find_library(MOSEK_LIB mosek REQUIRED)
if (MOSEK_INCLUDE_DIR AND MOSEK_LIB)
    message(STATUS "MOSEK found")
    set(MOSEK_LIBS ${MOSEK_LIB})
else()
    message(FATAL_ERROR "MOSEK library not found")
endif()

# FIND EIGEN
find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
    message(STATUS "Eigen found")
else()
    message(FATAL_ERROR "Eigen not found")
endif()

# FIND GSL
find_package(GSL REQUIRED)
if(GSL_FOUND)
    message(STATUS "GSL found")
else()
    message(FATAL_ERROR "GSL not found")
endif()

# Collect Sources
file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.c")

# Collect Header Directories
HEADER_DIRECTORIES(BFES_HEADER_DIRS)

# Executable
add_executable(BFES main_benchmarking.cpp ${SOURCES})

# Link Libraries
target_link_libraries(BFES OpenMP::OpenMP_CXX ${LIBOMP_LIBS} ${MOSEK_LIBS} ${GSL_LIBRARIES} Eigen3::Eigen)

# Include Headers
target_include_directories(BFES PUBLIC ${LIBOMP_INCLUDE_DIR} ${BFES_HEADER_DIRS} ${MOSEK_INCLUDE_DIR} ${GSL_INCLUDE_DIRS})