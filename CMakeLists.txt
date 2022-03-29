cmake_minimum_required(VERSION 3.17)
project(camel VERSION 0.2.0
              LANGUAGES C CXX
              DESCRIPTION "Snp detection tool")

# test if camel is main project
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(CAMEL_MAIN_PROJECT TRUE)
endif ()

# output dirs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# includes
include(cmake/Warnings.cmake)
include(GenerateExportHeader)
include(FetchContent)

# SKBUILD

# options
option(CAMEL_BUILD_SHARED_LIBS "Build camel and it's dependencies as shared lib" OFF)
option(CAMEL_BUILD_TESTS "Build tests" OFF)
option(CAMEL_PACK "Include camel packaging info" OFF)
option(CAMEL_BUILD_PYTHON "Build camel python bindings" OFF) 
option(CAMEL_BUILD_EXE "Build camel executable" OFF)
option(CAMEL_DEVEL_OPTS "Enable options used for development" OFF)

if (SKBUILD)
  set(CAMEL_BUILD_PYTHON ON)
endif ()

if (CAMEL_PACK AND SKBUILD)
  message(FATAL_ERROR "CAMEL_PACK and CAMEL_BUILD_PYTHON mutually exclusive")
  return()
endif ()

if (CAMEL_BUILD_PYTHON)
  set_property(GLOBAL # Sniff... Sniff... can do better?
               PROPERTY POSITION_INDEPENDENT_CODE TRUE)
endif ()

if (NOT DEFINED CMAKE_CXX_VISIBILITY_PRESET AND
    NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN AND NOT CAMEL_MAIN_PROJECT)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)
endif ()

# build type shared/static
if (CAMEL_BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON) 
endif ()

FetchContent_Declare(
  bioparser
  GIT_REPOSITORY https://github.com/rvaser/bioparser
  GIT_TAG 3.0.15) 

FetchContent_Declare(
  biosoup
  GIT_REPOSITORY https://github.com/rvaser/biosoup
  GIT_TAG 0.10.0)

FetchContent_Declare(
        cereal
        GIT_REPOSITORY https://github.com/USCiLab/cereal
        GIT_TAG v1.3.2)

FetchContent_GetProperties(cereal)
if (NOT cereal_POPULATED)
    FetchContent_Populate(cereal)
    add_subdirectory(
            ${cereal_SOURCE_DIR}
            ${cereal_BINARY_DIR}
            EXCLUDE_FROM_ALL)
    add_library(cereal::cereal ALIAS cereal)
endif ()

FetchContent_Declare(
  edlib
  GIT_REPOSITORY https://github.com/Martinsos/edlib
  GIT_TAG v1.2.7)

FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 8.1.1)

FetchContent_Declare(
  ram
  GIT_REPOSITORY https://github.com/lbcb-sci/ram
  GIT_TAG 2.1.1)

FetchContent_Declare(
  thread_pool
  GIT_REPOSITORY https://github.com/rvaser/thread_pool.git
  GIT_TAG 4.0.0)

FetchContent_MakeAvailable(bioparser biosoup cereal edlib fmt thread_pool ram)

# create main target
set(camel_SOURCES
  ${PROJECT_SOURCE_DIR}/src/camel/detail/overlap.cc

  ${PROJECT_SOURCE_DIR}/src/camel/coverage.cc
  ${PROJECT_SOURCE_DIR}/src/camel/mapping.cc
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
    biosoup thread_pool
  PRIVATE
    bioparser cereal edlib fmt ram)

# enable development options
if (CAMEL_DEVEL_OPTS)
  add_compile_options(
    $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>
    $<$<CONFIG:Debug,RelWithDebInfo>:-fno-omit-frame-pointer>)
  add_link_options(
    $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>)
endif ()

# generate export header and attach it to somelib
generate_export_header(camel EXPORT_FILE_NAME include/camel/export.h)
target_compile_definitions(camel PUBLIC
  "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:CAMEL_STATIC_DEFINE>")
target_include_directories(camel PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>")

if (CAMEL_BUILD_EXE)
  FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
    GIT_TAG v3.0.0)
    FetchContent_MakeAvailable(cxxopts)
 
  add_executable(camel_covg_exe ${PROJECT_SOURCE_DIR}/src/main.cc) 
  target_link_libraries(camel_covg_exe 
    PRIVATE 
      camel cxxopts fmt)
endif ()

if (CAMEL_BUILD_PYTHON)
  if (SKBUILD)
    message(STATUS "Python_VERSION ${PYTHON_VERSION_STRING}")
    message(STATUS "Python_EXECUTABLE ${PYTHON_EXECUTABLE}")
    message(STATUS "Python_INCLUDE_DIR ${PYTHON_INCLUDE_DIR}")
    message(STATUS "Python_LIBRARIES ${PYTHON_LIBRARY}")
    set(Python_VERSION "${PYTHON_VERSION_STRING}")
    set(Python_EXECUTABLE "${PYTHON_EXECUTABLE}")
    set(Python_INCLUDE_DIR "${PYTHON_INCLUDE_DIR}")
    set(Python_LIBRARIES "${PYTHON_LIBRARY}")
  endif ()

  find_package(Python COMPONENTS Interpreter Development REQUIRED)
  FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG master)

  FetchContent_MakeAvailable(nanobind)
  nanobind_add_module(camelpy_ext NB_STATIC 
    ${PROJECT_SOURCE_DIR}/src/bind.cc) 
  target_link_libraries(camelpy_ext PRIVATE camel thread_pool)

  if (SKBUILD)
    install(TARGETS camelpy_ext LIBRARY DESTINATION .)
  endif ()
elseif (CAMEL_PACK) # include packing for pure c++ install
  add_subdirectory(${PROJECT_SOURCE_DIR}/packaging)
endif ()

if (CAMEL_BUILD_TESTS)
  include(${PROJECT_SOURCE_DIR}/test/CamelTest.cmake)
endif ()
