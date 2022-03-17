cmake_minimum_required(VERSION 3.16)
project(camel VERSION 1.0.0
              LANGUAGES C CXX
              DESCRIPTION "Snp detection tool")

# includes
include(cmake/Warnings.cmake)
include(GenerateExportHeader)

# test if camel is main project
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(CAMEL_MAIN_PROJECT TRUE)
endif ()

# options
option(CAMEL_BUILD_SHARED_LIBS "Build camel and it's dependencies as shared lib" OFF)
option(CAMEL_BUILD_TESTS "Build tests" ${CAMEL_MAIN_PROJECT})
option(CAMEL_INCLUDE_PACKING "Inclue camel packaging info" ${CAMEL_MAIN_PROJECT})

if (NOT DEFINED CMAKE_CXX_VISIBILITY_PRESET AND
    NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)
endif ()

# build type shared/static
if (CAMEL_BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON) 
endif ()

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

# generate export header and attach it to somelib
generate_export_header(camel EXPORT_FILE_NAME include/camel/export.h)
target_compile_definitions(camel PUBLIC
  "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:CAMEL_STATIC_DEFINE>")
target_include_directories(camel PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>")

# include packaging
if (CAMEL_INCLUDE_PACKING)
  add_subdirectory(${PROJECT_SOURCE_DIR}/packaging)
endif ()
