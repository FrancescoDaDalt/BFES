cmake_minimum_required(VERSION 3.5.0)
project(BFES VERSION 0.1.0 LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


# Check if the compiler supports LTO
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-flto" COMPILER_SUPPORTS_LTO)

# Enable LTO if supported
if(COMPILER_SUPPORTS_LTO)
    # Add LTO flags for compiler
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto")
    
    # Add LTO flags for linker
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -flto")
else()
    message(WARNING "LTO is not supported by the compiler")
endif()

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

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -Xclang -fopenmp ")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -march=native -Xclang -fopenmp ")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(APPLE)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)

    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()
endif()

if (NOT MOSEK_INCLUDE_DIR)
    set(MOSEK_INCLUDE_DIR "/opt/mosek/10.0/tools/platform/osxaarch64/h")
endif()
if (NOT MOSEK_LIB)
    set(MOSEK_LIB "/opt/mosek/10.0/tools/platform/osxaarch64/bin/libfusion64.dylib")
endif()
# if (NOT MOSEK_LIB_DIR)
#     set(MOSEK_LIB_DIR "/opt/mosek/10.0/tools/platform/osxaarch64/bin")
# endif()
if (NOT LIBOMP_INCLUDE_DIR)
    set(LIBOMP_INCLUDE_DIR "/opt/homebrew/Cellar/libomp/19.1.5/include")
endif()
if (NOT LIBOMP_LIB)
    set(LIBOMP_LIB "/opt/homebrew/Cellar/libomp/19.1.5/lib/libomp.dylib") 
endif()

# set(CMAKE_INSTALL_RPATH "${MOSEK_LIB_DIR}")
# set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)


# FIND OPENMP
find_path(LIBOMP_INCLUDE_DIR omp.h REQUIRED)
find_library(LIBOMP_LIB libomp REQUIRED)
if (LIBOMP_INCLUDE_DIR AND LIBOMP_LIB)
    message(STATUS "LIBOMP found")
    include_directories(${LIBOMP_INCLUDE_DIR})
    set(LIBOMP_LIBS ${LIBOMP_LIB})
else()
    message(FATAL_ERROR "LIBOMP library not found")
endif()

# FIND MOSEK
find_path(MOSEK_INCLUDE_DIR mosek.h REQUIRED)
find_library(MOSEK_LIB mosek REQUIRED)
if (MOSEK_INCLUDE_DIR AND MOSEK_LIB)
    message(STATUS "MOSEK found")
    include_directories(${MOSEK_INCLUDE_DIR})
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

file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.c")
# file(GLOB_RECURSE HEADERS BFES/*.hpp)


link_directories(BFES ${MOSEK_LIB_DIR})

add_executable(BFES main_benchmarking.cpp ${SOURCES})

HEADER_DIRECTORIES(BFES_HEADER_DIRS)

# message(${BFES_HEADER_DIRS})
# message(${OpenMP_CXX_INCLUDE_DIRS})
# message(STATUS "OpenMP libraries: ${OpenMP_CXX_LIBRARIES}")
message(${MOSEK_LIBS})
message(${LIBOMP_LIBS})

target_link_libraries(BFES ${LIBOMP_LIBS} ${MOSEK_LIBS} ${GSL_LIBRARIES} Eigen3::Eigen)

target_include_directories(BFES PUBLIC ${BFES_HEADER_DIRS} ${OpenMP_CXX_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})
# target_include_directories(BFES PRIVATE BFES)