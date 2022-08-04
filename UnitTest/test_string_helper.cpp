#include <array>
#include <cstddef>
#include <vector>
#include "common/types.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "util/string_helper.h"

TEST(strStreamify, HandlesZeroInput) {
  using namespace std;
  using namespace testing;

  vector<GDouble> testDouble{0.0, 0.1234, 0.35632};
  array<GInt, 2>  testInt{1, -122334};

  ASSERT_EQ(strStreamify<3>(testDouble).str(), "0 0.1234 0.35632");
  ASSERT_EQ(strStreamify<2>(testInt).str(), "1 -122334");
}

TEST(toStringVector, HandlesZeroInput) {
  using namespace std;
  using namespace testing;

  vector<GDouble> testDouble{0.0, 0.1234, 0.35632};
  vector<byte>    testByte{static_cast<byte>(0), static_cast<byte>(1), static_cast<byte>(42)};

  ASSERT_THAT(toStringVector(testDouble), ElementsAre("0.000000", "0.123400", "0.356320"));
  ASSERT_THAT(toStringVector(testByte), ElementsAre("0", "1", "42"));
}