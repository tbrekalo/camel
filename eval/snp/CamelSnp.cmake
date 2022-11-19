find_package(nlohmann_json)

add_executable(snp_eval ${CMAKE_CURRENT_LIST_DIR}/src/main.cc)
target_link_libraries(
  snp_eval PRIVATE camel cxxopts::cxxopts edlib fmt::fmt
                   nlohmann_json::nlohmann_json tsl::robin_map)
