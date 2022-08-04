#include "gtest/gtest.h"
#include "math/hilbert.h"

TEST(HilbertIndex2D, HandlesZeroInput) {
  VectorD<2> quadrant0 = {0, 0};
  EXPECT_EQ(hilbert::index<2>(quadrant0, 1), 0);
  EXPECT_EQ(hilbert::index<2>(quadrant0, 2), 0);
  EXPECT_EQ(hilbert::index<2>(quadrant0, 3), 0);
  EXPECT_EQ(hilbert::index<2>(quadrant0, 4), 0);
  VectorD<2> quadrant1 = {0, 1};
  EXPECT_EQ(hilbert::index<2>(quadrant1, 1), 1);
  // level1: 1 * 4 level 2: 1
  EXPECT_EQ(hilbert::index<2>(quadrant1, 2), 5);
  EXPECT_EQ(hilbert::index<2>(quadrant1, 3), 21);
  EXPECT_EQ(hilbert::index<2>(quadrant1, 4), 85);
  VectorD<2> quadrant2 = {1, 1};
  EXPECT_EQ(hilbert::index<2>(quadrant2, 1), 2);
  // level1: 2 * 4 level 2: 2 -> 10
  EXPECT_EQ(hilbert::index<2>(quadrant2, 2), 10);
  EXPECT_EQ(hilbert::index<2>(quadrant2, 3), 42);
  EXPECT_EQ(hilbert::index<2>(quadrant2, 4), 170);
  VectorD<2> quadrant3 = {1, 0};
  EXPECT_EQ(hilbert::index<2>(quadrant3, 1), 3);
  EXPECT_EQ(hilbert::index<2>(quadrant3, 2), 15);
  EXPECT_EQ(hilbert::index<2>(quadrant3, 3), 63);
  EXPECT_EQ(hilbert::index<2>(quadrant3, 4), 255);

  // level1: 1 * 4 level 2: 2
  VectorD<2> quadrant6 = {0.3, 1};
  EXPECT_EQ(hilbert::index<2>(quadrant6, 2), 6);

  VectorD<2> quadrant5 = {0.2, 1};
  EXPECT_EQ(hilbert::index<2>(quadrant5, 2), 5);

  // level1: 2 * 4 level 2: 1 -> 9
  VectorD<2> quadrant9 = {0.5, 1};
  EXPECT_EQ(hilbert::index<2>(quadrant9, 2), 9);

  // level1: 1 * 4 level 2: 3 -> 7
  VectorD<2> quadrant7 = {0.4, 0.6};
  EXPECT_EQ(hilbert::index<2>(quadrant7, 2), 7);

  // level1: 2 * 4 level 2: 0 -> 8
  VectorD<2> quadrant8 = {0.5, 0.5};
  EXPECT_EQ(hilbert::index<2>(quadrant8, 2), 8);
}

TEST(HilbertIndex3D, HandlesZeroInput) {
  VectorD<3> quadrant0 = {0, 0, 0};
  EXPECT_EQ(hilbert::index<3>(quadrant0, 1), 0);
  EXPECT_EQ(hilbert::index<3>(quadrant0, 2), 0);
  EXPECT_EQ(hilbert::index<3>(quadrant0, 3), 0);
  EXPECT_EQ(hilbert::index<3>(quadrant0, 4), 0);
  VectorD<3> quadrant1 = {0, 1, 0};
  EXPECT_EQ(hilbert::index<3>(quadrant1, 1), 1);
  EXPECT_EQ(hilbert::index<3>(quadrant1, 2), 9);
  EXPECT_EQ(hilbert::index<3>(quadrant1, 3), 73);
  EXPECT_EQ(hilbert::index<3>(quadrant1, 4), 585);
  VectorD<3> quadrant2 = {1, 1, 0};
  EXPECT_EQ(hilbert::index<3>(quadrant2, 1), 2);
  EXPECT_EQ(hilbert::index<3>(quadrant2, 2), 18);
  EXPECT_EQ(hilbert::index<3>(quadrant2, 3), 146);
  EXPECT_EQ(hilbert::index<3>(quadrant2, 4), 1170);
  VectorD<3> quadrant3 = {1, 0, 0};
  EXPECT_EQ(hilbert::index<3>(quadrant3, 1), 3);
  EXPECT_EQ(hilbert::index<3>(quadrant3, 2), 27);
  EXPECT_EQ(hilbert::index<3>(quadrant3, 3), 219);
  EXPECT_EQ(hilbert::index<3>(quadrant3, 4), 1755);

  VectorD<3> quadrant4 = {1, 0, 1};
  EXPECT_EQ(hilbert::index<3>(quadrant4, 1), 4);
  EXPECT_EQ(hilbert::index<3>(quadrant4, 2), 36);
  EXPECT_EQ(hilbert::index<3>(quadrant4, 3), 292);
  EXPECT_EQ(hilbert::index<3>(quadrant4, 4), 2340);

  VectorD<3> quadrant5 = {0, 0, 1};
  EXPECT_EQ(hilbert::index<3>(quadrant5, 1), 5);
  EXPECT_EQ(hilbert::index<3>(quadrant5, 2), 45);
  EXPECT_EQ(hilbert::index<3>(quadrant5, 3), 365);
  EXPECT_EQ(hilbert::index<3>(quadrant5, 4), 2925);

  VectorD<3> quadrant6 = {0, 1, 1};
  EXPECT_EQ(hilbert::index<3>(quadrant6, 1), 6);
  EXPECT_EQ(hilbert::index<3>(quadrant6, 2), 54);
  EXPECT_EQ(hilbert::index<3>(quadrant6, 3), 438);
  EXPECT_EQ(hilbert::index<3>(quadrant6, 4), 3510);

  VectorD<3> quadrant7 = {1, 1, 1};
  EXPECT_EQ(hilbert::index<3>(quadrant7, 1), 7);
  EXPECT_EQ(hilbert::index<3>(quadrant7, 2), 63);
  EXPECT_EQ(hilbert::index<3>(quadrant7, 3), 511);
  EXPECT_EQ(hilbert::index<3>(quadrant7, 4), 4095);

  // quadrants with close coordinates!
  //  quadrant spacing level 1 0.5 <- identical id on level 1
  //  quadrant spacing level 2 0.25 <- different id on level 2
  VectorD<3> quadrantA = {0.375, 0.125, 0.625};
  VectorD<3> quadrantB = {0.375, 0.125, 0.875};
  EXPECT_EQ(hilbert::index<3>(quadrantA, 1), 5);
  EXPECT_EQ(hilbert::index<3>(quadrantA, 2), 43);
  EXPECT_EQ(hilbert::index<3>(quadrantB, 1), 5);
  EXPECT_EQ(hilbert::index<3>(quadrantB, 2), 44);

  VectorD<3> quadrantC = {0.875, 0.125, 0.875};
  VectorD<3> quadrantD = {0.625, 0.375, 0.875};
  EXPECT_EQ(hilbert::index<3>(quadrantC, 1), 4);
  EXPECT_EQ(hilbert::index<3>(quadrantC, 2), 36);
  EXPECT_EQ(hilbert::index<3>(quadrantD, 1), 4);
  EXPECT_EQ(hilbert::index<3>(quadrantD, 2), 38);
}

TEST(HilbertIndex4D, HandlesZeroInput) {
  VectorD<4> quadrant0 = {0, 0, 0, 0};
  EXPECT_EQ(hilbert::index<4>(quadrant0, 1), 0);
  EXPECT_EQ(hilbert::index<4>(quadrant0, 2), 0);
  EXPECT_EQ(hilbert::index<4>(quadrant0, 3), 0);
  EXPECT_EQ(hilbert::index<4>(quadrant0, 4), 0);

  VectorD<4> quadrant1 = {0, 1, 0, 0};
  EXPECT_EQ(hilbert::index<4>(quadrant1, 1), 1);
  EXPECT_EQ(hilbert::index<4>(quadrant1, 2), 17);
  EXPECT_EQ(hilbert::index<4>(quadrant1, 3), 273);
  EXPECT_EQ(hilbert::index<4>(quadrant1, 4), 4369);

  VectorD<4> quadrant2 = {1, 1, 0, 0};
  EXPECT_EQ(hilbert::index<4>(quadrant2, 1), 2);
  EXPECT_EQ(hilbert::index<4>(quadrant2, 2), 34);
  EXPECT_EQ(hilbert::index<4>(quadrant2, 3), 546);
  EXPECT_EQ(hilbert::index<4>(quadrant2, 4), 8738);

  VectorD<4> quadrant3 = {1, 0, 0, 0};
  EXPECT_EQ(hilbert::index<4>(quadrant3, 1), 3);
  EXPECT_EQ(hilbert::index<4>(quadrant3, 2), 51);
  EXPECT_EQ(hilbert::index<4>(quadrant3, 3), 819);
  EXPECT_EQ(hilbert::index<4>(quadrant3, 4), 13107);

  VectorD<4> quadrant4 = {1, 0, 1, 0};
  EXPECT_EQ(hilbert::index<4>(quadrant4, 1), 4);
  EXPECT_EQ(hilbert::index<4>(quadrant4, 2), 68);
  EXPECT_EQ(hilbert::index<4>(quadrant4, 3), 1092);
  EXPECT_EQ(hilbert::index<4>(quadrant4, 4), 17476);

  VectorD<4> quadrant5 = {0, 0, 1, 0};
  EXPECT_EQ(hilbert::index<4>(quadrant5, 1), 5);
  EXPECT_EQ(hilbert::index<4>(quadrant5, 2), 85);
  EXPECT_EQ(hilbert::index<4>(quadrant5, 3), 1365);
  EXPECT_EQ(hilbert::index<4>(quadrant5, 4), 21845);

  VectorD<4> quadrant6 = {0, 1, 1, 0};
  EXPECT_EQ(hilbert::index<4>(quadrant6, 1), 6);
  EXPECT_EQ(hilbert::index<4>(quadrant6, 2), 102);
  EXPECT_EQ(hilbert::index<4>(quadrant6, 3), 1638);
  EXPECT_EQ(hilbert::index<4>(quadrant6, 4), 26214);

  VectorD<4> quadrant7 = {1, 1, 1, 0};
  EXPECT_EQ(hilbert::index<4>(quadrant7, 1), 7);
  EXPECT_EQ(hilbert::index<4>(quadrant7, 2), 119);
  EXPECT_EQ(hilbert::index<4>(quadrant7, 3), 1911);
  EXPECT_EQ(hilbert::index<4>(quadrant7, 4), 30583);

  // quadrants with close coordinates!
  //  quadrant spacing level 1 0.5 <- identical id on level 1
  //  quadrant spacing level 2 0.25 <- different id on level 2
  VectorD<4> quadrantA = {0.375, 0.125, 0.625, 0.875};
  VectorD<4> quadrantB = {0.375, 0.125, 0.875, 0.875};

  EXPECT_EQ(hilbert::index<4>(quadrantA, 1), 15);
  EXPECT_EQ(hilbert::index<4>(quadrantA, 2), 249);
  EXPECT_EQ(hilbert::index<4>(quadrantB, 1), 15);
  EXPECT_EQ(hilbert::index<4>(quadrantB, 2), 254);

  VectorD<4> quadrantC = {0.875, 0.125, 0.875, 0.375};
  VectorD<4> quadrantD = {0.625, 0.375, 0.875, 0.375};
  VectorD<4> quadrantE = {0.125, 0.625, 0.625, 0.125};
  VectorD<4> quadrantF = {0.625, 0.125, 0.875, 0.375};

  EXPECT_EQ(hilbert::index<4>(quadrantD, 1), 4);
  EXPECT_EQ(hilbert::index<4>(quadrantD, 2), 76);
  EXPECT_EQ(hilbert::index<4>(quadrantE, 1), 6);
  EXPECT_EQ(hilbert::index<4>(quadrantE, 2), 96);
  EXPECT_EQ(hilbert::index<4>(quadrantC, 1), 4);
  EXPECT_EQ(hilbert::index<4>(quadrantC, 2), 78);
  EXPECT_EQ(hilbert::index<4>(quadrantF, 1), 4);
  EXPECT_EQ(hilbert::index<4>(quadrantF, 2), 79);

  VectorD<4> quadrantG = {0.625, 0.125, 0.875, 0.875};
  VectorD<4> quadrantH = {0.875, 0.125, 0.875, 0.875};
  EXPECT_EQ(hilbert::index<4>(quadrantG, 1), 14);
  EXPECT_EQ(hilbert::index<4>(quadrantG, 2), 239);
  EXPECT_EQ(hilbert::index<4>(quadrantH, 1), 14);
  EXPECT_EQ(hilbert::index<4>(quadrantH, 2), 238);

  VectorD<4> quadrantI = {0.625, 0.625, 0.875, 0.875};
  VectorD<4> quadrantJ = {0.625, 0.875, 0.875, 0.875};
  EXPECT_EQ(hilbert::index<4>(quadrantI, 1), 13);
  EXPECT_EQ(hilbert::index<4>(quadrantI, 2), 223);
  EXPECT_EQ(hilbert::index<4>(quadrantJ, 1), 13);
  EXPECT_EQ(hilbert::index<4>(quadrantJ, 2), 220);

  VectorD<4> quadrantK = {0.625, 0.625, 0.875, 0.375};
  VectorD<4> quadrantL = {0.875, 0.625, 0.875, 0.375};
  EXPECT_EQ(hilbert::index<4>(quadrantK, 1), 7);
  EXPECT_EQ(hilbert::index<4>(quadrantK, 2), 127);
  EXPECT_EQ(hilbert::index<4>(quadrantL, 1), 7);
  EXPECT_EQ(hilbert::index<4>(quadrantL, 2), 126);

  VectorD<4> quadrantM = {0.125, 0.625, 0.875, 0.375};
  VectorD<4> quadrantN = {0.375, 0.625, 0.875, 0.375};
  EXPECT_EQ(hilbert::index<4>(quadrantM, 1), 6);
  EXPECT_EQ(hilbert::index<4>(quadrantM, 2), 111);
  EXPECT_EQ(hilbert::index<4>(quadrantN, 1), 6);
  EXPECT_EQ(hilbert::index<4>(quadrantN, 2), 110);

  VectorD<4> quadrantO = {0.125, 0.125, 0.875, 0.375};
  VectorD<4> quadrantP = {0.375, 0.125, 0.875, 0.375};
  EXPECT_EQ(hilbert::index<4>(quadrantO, 1), 5);
  EXPECT_EQ(hilbert::index<4>(quadrantO, 2), 95);
  EXPECT_EQ(hilbert::index<4>(quadrantP, 1), 5);
  EXPECT_EQ(hilbert::index<4>(quadrantP, 2), 94);
}
