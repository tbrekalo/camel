include(FetchContent)

find_package(Catch2 3 QUIET)
if (NOT Catch2_FOUND)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.0.0-preview4
  )

  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()

set(camel_TESTS_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/src/serialization.cc)

add_executable(camel_tests ${${PROJECT_NAME}_TESTS_SOURCES})
target_link_libraries(camel_tests
  PRIVATE
    camel
    Catch2::Catch2WithMain)

include(CTest)
include(Catch)
enable_testing()
catch_discover_tests(camel_tests
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test)
