#include "camel/coverage.h"
#include "camel/io.h"
#include "camel/state.h"
#include "catch2/catch_test_macros.hpp"

// TODO: snif... snif...
std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects = 0U;

namespace camel::test {

static auto const kTestPile00 = Pile{0U, "test_pile_00",
                                     std::vector<Coverage>{
                                         Coverage{1, 0, 0, 0, 0, 0},
                                         Coverage{0, 1, 0, 0, 0, 0},
                                         Coverage{0, 0, 1, 0, 0, 0},
                                         Coverage{0, 0, 0, 1, 0, 0},
                                         Coverage{0, 0, 0, 0, 1, 0},
                                         Coverage{0, 0, 0, 0, 0, 1},
                                     }};

static auto const kTestPile01 = Pile{0U, "test_pile_01",
                                     std::vector<Coverage>{
                                         Coverage{0, 0, 0, 0, 0, 1},
                                         Coverage{0, 0, 0, 0, 1, 0},
                                         Coverage{0, 0, 0, 1, 0, 0},
                                         Coverage{0, 0, 1, 0, 0, 0},
                                         Coverage{0, 1, 0, 0, 0, 0},
                                         Coverage{1, 0, 0, 0, 0, 0},
                                     }};

static auto const kTestDumpPath = std::filesystem::path("./test_dump");

}  // namespace camel::test

TEST_CASE("camel pile serialization", "[camel][pile][coverage][serialize]") {
  auto state =
      camel::State{.thread_pool = std::make_shared<thread_pool::ThreadPool>(1),
                   .log_path = "./temp_camel_test_log"};

  if (std::filesystem::exists(state.log_path)) {
    std::filesystem::remove(state.log_path);
  }
  std::filesystem::create_directory(state.log_path);

  auto const src_piles = std::vector<camel::Pile>{camel::test::kTestPile00,
                                                  camel::test::kTestPile01};

  camel::SerializePiles(state, src_piles, camel::test::kTestDumpPath);

  auto const read_piles =
      camel::DeserializePiles(state, camel::test::kTestDumpPath);

  REQUIRE(src_piles.size() == read_piles.size());
  for (auto i = 0UL; i < src_piles.size(); ++i) {
    CHECK(src_piles[i].id == read_piles[i].id);
    CHECK(src_piles[i].seq_name == src_piles[i].seq_name);

    REQUIRE(src_piles[i].covgs.size() == read_piles[i].covgs.size());
    for (auto j = 0UL; j < src_piles.size(); ++j) {
      CHECK(src_piles[i].covgs[j].a == read_piles[i].covgs[j].a);
      CHECK(src_piles[i].covgs[j].c == read_piles[i].covgs[j].c);
      CHECK(src_piles[i].covgs[j].g == read_piles[i].covgs[j].g);
      CHECK(src_piles[i].covgs[j].t == read_piles[i].covgs[j].t);
      CHECK(src_piles[i].covgs[j].del == read_piles[i].covgs[j].del);
      CHECK(src_piles[i].covgs[j].ins == read_piles[i].covgs[j].ins);
    }
  }
}
