include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if (NOT DEFINED camel_INSTALL_CMAKEDIR) 
  set(camel_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/camel"
    CACHE STRING "Path to camel cmake files")
endif ()

install(TARGETS camel EXPORT camel_TARGETS
        RUNTIME COMPONENT camel_Runtime
        LIBRARY COMPONENT camel_Runtime
        NAMELINK_COMPONENT camel_Development
        ARCHIVE COMPONENT camel_Development
        INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(DIRECTORY "${camel_SOURCE_DIR}/include/" "${camel_BINARY_DIR}/include/"
        TYPE INCLUDE
        COMPONENT camel_Development)

if (BUILD_SHARED_LIBS)
  set(type shared)
else ()
  set(type static)
endif ()

install(EXPORT camel_TARGETS
        DESTINATION "${camel_INSTALL_CMAKEDIR}"
        NAMESPACE camel
        FILE camel-${type}-targets.cmake
        COMPONENT camel_Development)

write_basic_package_version_file(
  camelConfigVersion.cmake 
  COMPATIBILITY SameMajorVersion)

install(FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/camelConfigVersion.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/camelConfigVersion.cmake"
        DESTINATION "${camel_INSTALL_CMAKEDIR}"
        COMPONENT camel_Development)

include(CPack)
