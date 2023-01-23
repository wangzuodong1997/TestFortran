MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

set (MKLROOT "C:/Program Files (x86)/Intel/oneAPI/mkl/latest")
message(STATUS "MKLROOT: " ${MKLROOT})

# Enable Fortran
set(CMAKE_Fortran_COMPILER "C:/Program Files (x86)/Intel/oneAPI/compiler/2023.0.0/windows/bin/intel64/ifort.exe")
enable_language(Fortran)


# Project name and version
project(stokes)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008")
    set(bounds "-Wall -fbounds-check -fimplicit-none -C -fcheck=all ")
  endif()

#if (Fortran_COMPILER_NAME MATCHES "f95")
#  set(bounds "-Wall -fimplicit-none -fbounds-check")
#endif()
  

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_RELEASE " -O3 " CACHE STRING
  " Flags used by the compiler during release builds. ")
set(CMAKE_Fortran_FLAGS_NATIVE " -O4 -march=native -mtune=native -Ofast " CACHE STRING
  " Flags used by the compiler during release builds. ")

# Define runtime output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../EXECUTABLE)

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

ADD_CUSTOM_TARGET(native
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Native ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Native"
  )

# Since the original sources directory *should* not change structure, i.e.
# files are not added/removed, we can use the file()
# command to populate a variable with the list of needed files

file(GLOB_RECURSE global_sources ${CMAKE_SOURCE_DIR}/SOURCES/*.f90)
#message(STATUS "cmakedir: ${CMAKE_SOURCE_DIR}")
#Remove items



# add executable
add_executable(a.exe
  ${global_sources}
  )

# find Intel MKL
#set(BLA_VENDOR Intel10_64lp_seq)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)



# include and link library

target_include_directories(a.exe PUBLIC ${MKLROOT}/include)
target_link_libraries(a.exe LAPACK::LAPACK BLAS::BLAS )


