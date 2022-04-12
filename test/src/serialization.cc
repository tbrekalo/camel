#include "camel/coverage.h"
#include "camel/io.h"
#include "catch2/catch_test_macros.hpp"

// TODO: snif... snif...
std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects = 0U;

namespace camel::test {

static auto const kTestPile00 = Pile(0U, "test_pile_00",
                                     std::vector<Coverage>{
                                         Coverage(1, 0, 0, 0),
                                         Coverage(0, 1, 0, 0),
                                         Coverage(0, 0, 1, 0),
                                         Coverage(0, 0, 0, 1),
                                     });

static auto const kTestPile01 = Pile{0U, "test_pile_01",
                                     std::vector<Coverage>{
                                         Coverage(0, 0, 0, 1),
                                         Coverage(0, 0, 1, 0),
                                         Coverage(0, 1, 0, 0),
                                         Coverage(1, 0, 0, 0),
                                     }};

static auto const kTestDumpPath = std::filesystem::path("./test_dump");

}  // namespace camel::test

TEST_CASE("camel pile serialization", "[camel][pile][coverage][serialize]") {
  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(1);
  auto const src_piles = std::vector<camel::Pile>{camel::test::kTestPile00,
                                                  camel::test::kTestPile01};

  camel::SerializePiles(thread_pool, src_piles, camel::test::kTestDumpPath);

  auto const read_piles =
      camel::DeserializePiles(thread_pool, camel::test::kTestDumpPath);

  REQUIRE(src_piles.size() == read_piles.size());
  for (auto i = 0UL; i < src_piles.size(); ++i) {
    CHECK(src_piles[i].id == read_piles[i].id);
    CHECK(src_piles[i].seq_name == src_piles[i].seq_name);

    REQUIRE(src_piles[i].covgs.size() == read_piles[i].covgs.size());
    for (auto j = 0UL; j < src_piles.size(); ++j) {
      CHECK(src_piles[i].covgs[j].mat == read_piles[i].covgs[j].mat);
      CHECK(src_piles[i].covgs[j].del == read_piles[i].covgs[j].del);
      CHECK(src_piles[i].covgs[j].ins == read_piles[i].covgs[j].ins);
      CHECK(src_piles[i].covgs[j].mis == read_piles[i].covgs[j].mis);
    }
  }
}
