cmake_minimum_required(VERSION 3.17)
project(camel VERSION 0.8.0
              LANGUAGES CXX
              DESCRIPTION "Read polishing tool")

# test if camel is main project
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(CAMEL_MAIN_PROJECT TRUE)
endif ()

# output dirs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# get conan
if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
  message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
  file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/0.18.1/conan.cmake"
                "${CMAKE_BINARY_DIR}/conan.cmake"
                TLS_VERIFY ON)
endif()

# includes
include(cmake/Warnings.cmake)
include(GenerateExportHeader)
include(FetchContent)

include(${CMAKE_BINARY_DIR}/conan.cmake)

# prepend conan path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})

# options
option(CAMEL_BUILD_SHARED_LIBS "Build camel and it's dependencies as shared lib" OFF)
option(CAMEL_BUILD_TESTS "Build tests" OFF)
option(CAMEL_PACK "Include camel packaging info" OFF)
option(CAMEL_BUILD_EXE "Build camel executable" OFF)
option(CAMEL_DEVEL_OPTS "Enable options used for development" OFF)

if (NOT DEFINED CMAKE_CXX_VISIBILITY_PRESET AND
    NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN AND NOT CAMEL_MAIN_PROJECT)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)
endif ()

# build type shared/static
if (CAMEL_BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON) 
endif ()

conan_cmake_configure(
  REQUIRES 
    cxxopts/3.0.0 
    fmt/8.1.1
    onetbb/2021.3.0
    span-lite/0.10.3
    tsl-robin-map/1.0.1
  GENERATORS 
    cmake_find_package)

conan_cmake_autodetect(settings)
conan_cmake_install(PATH_OR_REFERENCE .
                    BUILD missing
                    REMOTE conancenter
                    SETTINGS ${settings})

find_package(fmt)
find_package(span-lite)
find_package(tsl-robin-map)
find_package(TBB)

FetchContent_Declare(
  bioparser
  GIT_REPOSITORY https://github.com/rvaser/bioparser
  GIT_TAG 3.0.15) 

FetchContent_Declare(
  biosoup
  GIT_REPOSITORY https://github.com/rvaser/biosoup
  GIT_TAG 0.10.0)

FetchContent_Declare(
  edlib
  GIT_REPOSITORY https://github.com/Martinsos/edlib
  GIT_TAG v1.2.7)

FetchContent_Declare(
  spoa
  GIT_REPOSITORY https://github.com/rvaser/spoa
  GIT_TAG 4.0.7)

FetchContent_MakeAvailable(bioparser biosoup edlib spoa)

# create main target
set(camel_SOURCES
  ${PROJECT_SOURCE_DIR}/src/camel/detail/alignment.cc
  ${PROJECT_SOURCE_DIR}/src/camel/detail/call_sites.cc
  ${PROJECT_SOURCE_DIR}/src/camel/detail/correct.cc
  ${PROJECT_SOURCE_DIR}/src/camel/detail/coverage.cc
  ${PROJECT_SOURCE_DIR}/src/camel/detail/interval.cc
  ${PROJECT_SOURCE_DIR}/src/camel/detail/nucleic_view.cc
  ${PROJECT_SOURCE_DIR}/src/camel/detail/overlap.cc
  ${PROJECT_SOURCE_DIR}/src/camel/detail/task_queue.cc
  ${PROJECT_SOURCE_DIR}/src/camel/detail/window.cc

  ${PROJECT_SOURCE_DIR}/src/camel/correct.cc
  ${PROJECT_SOURCE_DIR}/src/camel/io.cc)

add_library(camel ${camel_SOURCES})
add_library(${PROJECT_NAME}::camel ALIAS camel)
set_target_properties(camel PROPERTIES
                      VERSION ${camel_VERSION}
                      SOVERSION ${camel_VERSION_MAJOR})
target_include_directories(camel 
  PUBLIC 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_compile_features(camel PUBLIC cxx_std_17)
target_link_libraries(camel
  PUBLIC
    biosoup TBB::tbb
  PRIVATE
    bioparser edlib fmt::fmt nonstd::span-lite spoa tsl::robin_map)

# enable development options
if (CAMEL_DEVEL_OPTS)
  target_compile_options(camel PUBLIC
    $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>
    $<$<CONFIG:Debug,RelWithDebInfo>:-fno-omit-frame-pointer>)
  target_link_options(camel PUBLIC
    $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>)
endif ()

# generate export header and attach it to somelib
generate_export_header(camel EXPORT_FILE_NAME include/camel/export.h)
target_compile_definitions(camel PUBLIC
  "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:CAMEL_STATIC_DEFINE>")
target_include_directories(camel PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>")

if (CAMEL_BUILD_EXE) 
  find_package(cxxopts)

  add_executable(camel_exe ${PROJECT_SOURCE_DIR}/src/main.cc) 
  target_link_libraries(camel_exe
    PRIVATE 
      camel cxxopts::cxxopts 
      fmt::fmt)
endif ()
