cmake_minimum_required(VERSION 3.16)
project(camel VERSION 0.1.0
              LANGUAGES C CXX
              DESCRIPTION "Snp detection tool")

# includes
include(cmake/Warnings.cmake)
include(GenerateExportHeader)
include(FetchContent)

# test if camel is main project
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(CAMEL_MAIN_PROJECT TRUE)
endif ()

# options
option(CAMEL_BUILD_SHARED_LIBS "Build camel and it's dependencies as shared lib" OFF)
option(CAMEL_BUILD_TESTS "Build tests" ${CAMEL_MAIN_PROJECT})
option(CAMEL_INCLUDE_PACKING "Inclue camel packaging info" ${CAMEL_MAIN_PROJECT})
option(CAMEL_BUILD_EXE "Build camel executable" OFF)

if (NOT DEFINED CMAKE_CXX_VISIBILITY_PRESET AND
    NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN AND NOT CAMEL_MAIN_PROJECT)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)
endif ()

# build type shared/static
if (CAMEL_BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON) 
endif ()

# lib dependencies
find_package(biosoup 0.10.0 QUIET)
if(NOT biosoup_FOUND)
  FetchContent_Declare(
    biosoup
    GIT_REPOSITORY https://github.com/rvaser/biosoup
    GIT_TAG 0.10.0)

  FetchContent_MakeAvailable(biosoup)
endif ()

find_package(edlib 1.2.7 QUIET)
if (NOT edlib_FOUND)
  FetchContent_Declare(
    edlib
    GIT_REPOSITORY https://github.com/Martinsos/edlib
    GIT_TAG v1.2.7)

  FetchContent_MakeAvailable(edlib)
endif()

find_package(fmt 8 QUIET)
if (NOT fmt_FOUND)
  FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 8.1.1)
endif ()

find_package(ram 2.1.1 QUIET)
if (NOT ram_FOUND)
  FetchContent_Declare(
    ram
    GIT_REPOSITORY https://github.com/lbcb-sci/ram
    GIT_TAG 2.1.1)

  FetchContent_MakeAvailable(ram)
endif()

find_package(thread_pool 4 QUIET)
if (NOT thread_pool_FOUND)
  FetchContent_Declare(
    thread_pool
    GIT_REPOSITORY https://github.com/rvaser/thread_pool.git
    GIT_TAG v4.0.0)
endif ()


FetchContent_MakeAvailable(biosoup edlib fmt thread_pool ram)

# create main target
set(camel_SOURCES
  ${PROJECT_SOURCE_DIR}/src/camel.cc)

add_library(camel ${camel_SOURCES})
add_library(${PROJECT_NAME}::camel ALIAS camel)
set_target_properties(camel PROPERTIES
                      VERSION ${camel_VERSION}
                      SOVERSION ${camel_VERSION_MAJOR})
target_include_directories(camel PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_compile_features(camel PUBLIC cxx_std_17)
target_link_libraries(camel
  PUBLIC
    biosoup thread_pool
  PRIVATE
    edlib fmt ram)

# generate export header and attach it to somelib
generate_export_header(camel EXPORT_FILE_NAME include/camel/export.h)
target_compile_definitions(camel PUBLIC
  "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:CAMEL_STATIC_DEFINE>")
target_include_directories(camel PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>")

if (CAMEL_BUILD_EXE)
  find_package(bioparser 3.0.15 QUIET)
  if (NOT bioparser_FOUND)
    FetchContent_Declare(
      bioparser
      GIT_REPOSITORY https://github.com/rvaser/bioparser
      GIT_TAG 3.0.15
    ) 
  
    FetchContent_MakeAvailable(bioparser)
  endif()

  find_package(cxxopts 3 QUIET)
  if (NOT cxxopts_FOUND)
    FetchContent_Declare(
      cxxopts
      GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
      GIT_TAG v3.0.0)
  endif ()
  
  FetchContent_MakeAvailable(bioparser cxxopts)

  add_executable(camel_exe ${PROJECT_SOURCE_DIR}/src/main.cc) 
  target_link_libraries(camel_exe 
    PRIVATE 
      bioparser camel cxxopts fmt)
endif ()

# include packaging
if (CAMEL_INCLUDE_PACKING)
  add_subdirectory(${PROJECT_SOURCE_DIR}/packaging)
endif ()

if (CAMEL_BUILD_TESTS)
  include(${PROJECT_SOURCE_DIR}/test/CamelTest.cmake)
endif ()
