#include "camel/serialization.h"

#include "camel/coverage.h"
#include "camel/io.h"
#include "catch2/catch_test_macros.hpp"

// TODO: snif... snif...
std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects = 0U;

namespace camel::test {

static auto const kTestDumpPath = std::filesystem::path("./test_dump");

struct Blob {
  std::uint32_t val_a;
  std::uint8_t val_b;

  auto operator==(Blob const& that) const -> bool {
    return val_a == that.val_a && val_b == that.val_b;
  }
};

struct NamedBlobs {
  std::uint32_t id;
  std::string name;
  std::vector<Blob> blobs;

  auto operator==(NamedBlobs const& that) const -> bool {
    return id == that.id && name == that.name && blobs == that.blobs;
  }
};

template <class Buff>
auto Store(Buff& buff, NamedBlobs const& named_blobs) -> void {
  buff(named_blobs.id, named_blobs.name, named_blobs.blobs);
}

template <class Buff>
auto Load(Buff& buff, NamedBlobs& named_blobs) -> void {
  buff(named_blobs.id, named_blobs.name, named_blobs.blobs);
}

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

}  // namespace camel::test

TEST_CASE("Binary serialization trivial", "[serialization][binary][trivial]") {
  auto binary_out = camel::BinaryOutBuffer();
  auto const test_src_blob = camel::test::Blob{.val_a = 0, .val_b = 1};

  binary_out(test_src_blob);

  SECTION("convert binary_out to binary_in") {
    auto binary_in = camel::BinaryInBuffer(binary_out);
    auto test_dst_blob = camel::test::Blob{.val_a = 1, .val_b = 0};

    binary_in(test_dst_blob);

    REQUIRE(test_src_blob == test_dst_blob);
  }

  SECTION("serialize to file with zlib") {
    camel::GzStoreBytes(binary_out.Bytes(),
                        camel::test::kTestDumpPath / "test_blob.gzip");
    auto binary_in = camel::BinaryInBuffer(
        camel::GzLoadBytes(camel::test::kTestDumpPath / "test_blob.gzip"));
    auto test_dst_blob = camel::test::Blob{.val_a = 1, .val_b = 0};

    binary_in(test_dst_blob);

    REQUIRE(test_src_blob == test_dst_blob);
  }
}

TEST_CASE("Binary serialization trivial vector",
          "[serialization][binary][vector][trivial]") {
  auto binary_out = camel::BinaryOutBuffer();
  auto test_src_blobs = std::vector<camel::test::Blob>();
  std::generate_n(std::back_inserter(test_src_blobs), 10U,
                  [i = 0U]() mutable -> camel::test::Blob {
                    return camel::test::Blob{
                        .val_a = i,
                        .val_b = static_cast<std::uint8_t>(10U - ++i)};
                  });

  binary_out(test_src_blobs);

  auto binary_in = camel::BinaryInBuffer(binary_out);
  auto test_dst_blobs = std::vector<camel::test::Blob>();

  binary_in(test_dst_blobs);

  REQUIRE(test_src_blobs.size() == test_dst_blobs.size());
  for (auto i = 0; i < test_src_blobs.size(); ++i) {
    CHECK(test_dst_blobs[i] == test_dst_blobs[i]);
  }
}

TEST_CASE("Binary serialization string",
          "[serialization][binary][string][trivial]") {
  using namespace std::literals;

  auto binary_out = camel::BinaryOutBuffer();
  auto test_src_str = "Mahuna Matata\0Maybe?"s;

  binary_out(test_src_str);

  auto binary_in = camel::BinaryInBuffer(binary_out);
  auto test_dst_str = ""s;

  binary_in(test_dst_str);
  REQUIRE(test_src_str.size() == test_dst_str.size());
  REQUIRE(test_src_str == test_dst_str);
}

TEST_CASE("Binary serialization nontrivial",
          "[serialization][binary][nontrivial]") {
  using namespace std::literals;

  auto binary_out = camel::BinaryOutBuffer();
  auto test_src_nblob = camel::test::NamedBlobs{
      .id = 42,
      .name = "named_blob",
      .blobs = std::vector<camel::test::Blob>{
          camel::test::Blob{.val_a = 314, .val_b = 42U}}};

  binary_out(test_src_nblob);

  SECTION("convert binary_out to binary_in") {
    auto binary_in = camel::BinaryInBuffer(binary_out);
    auto test_dst_nblob = camel::test::NamedBlobs{};

    binary_in(test_dst_nblob);
    REQUIRE(test_src_nblob == test_dst_nblob);
  }

  SECTION("serialize to file with zlib") {
    camel::GzStoreBytes(binary_out.Bytes(),
                        camel::test::kTestDumpPath / "test_nblob.gzip");

    auto binary_in = camel::BinaryInBuffer(
        camel::GzLoadBytes(camel::test::kTestDumpPath / "test_nblob.gzip"));
    auto test_dst_nblob = camel::test::NamedBlobs();
    binary_in(test_dst_nblob);

    REQUIRE(test_src_nblob == test_dst_nblob);
  }
}

TEST_CASE("Binary serialization nontrivial vector",
          "[serialization][binary][string][vector][nontrivial]") {
  using namespace std::literals;

  auto binary_out = camel::BinaryOutBuffer();
  auto test_src_nblobs = std::vector<camel::test::NamedBlobs>();

  auto const gen_vec_of_blobs =
      [](std::uint32_t const first_id,
         std::uint32_t const last_id) -> std::vector<camel::test::Blob> {
    auto dst = std::vector<camel::test::Blob>();
    auto const sz = last_id - first_id;

    dst.reserve(sz);
    std::generate_n(std::back_inserter(dst), sz,
                    [curr_id = first_id]() mutable -> camel::test::Blob {
                      return camel::test::Blob{
                          .val_a = curr_id++,
                          .val_b = static_cast<std::uint8_t>(curr_id % 5U)};
                    });

    return dst;
  };

  std::generate_n(std::back_inserter(test_src_nblobs), 1U,
                  [&gen_vec_of_blobs,
                   curr_blob_id = 0U]() mutable -> camel::test::NamedBlobs {
                    auto const dst = camel::test::NamedBlobs{
                        .id = curr_blob_id,
                        .name = "dummy_name",
                        .blobs = gen_vec_of_blobs(curr_blob_id * 5,
                                                  (curr_blob_id + 1U) * 5)};
                    ++curr_blob_id;
                    return dst;
                  });

  binary_out(test_src_nblobs);

  SECTION("convert binary_out to binary_in") {
    auto binary_in = camel::BinaryInBuffer(binary_out);
    auto test_dst_nblobs = std::vector<camel::test::NamedBlobs>();

    binary_in(test_dst_nblobs);

    REQUIRE(test_src_nblobs.size() == test_dst_nblobs.size());
    REQUIRE(test_src_nblobs == test_dst_nblobs);
  }

  SECTION("serialize to file with zlib") {
    camel::GzStoreBytes(binary_out.Bytes(),
                        camel::test::kTestDumpPath / "test_nblobs.gzip");
    auto binary_in = camel::BinaryInBuffer(
        camel::GzLoadBytes(camel::test::kTestDumpPath / "test_nblobs.gzip"));

    auto test_dst_nblobs = std::vector<camel::test::NamedBlobs>();

    binary_in(test_dst_nblobs);

    REQUIRE(test_src_nblobs.size() == test_dst_nblobs.size());
    REQUIRE(test_src_nblobs == test_dst_nblobs);
  }
}

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
