#pragma once

#include <gtest/gtest.h>
#include <cmath>

#define ABS_TEST(a, b, tol) EXPECT_LE(std::abs(((a) - (b))), (tol))

#define REL_TEST(a, b, tol) \
  if (std::abs(b) < 1e-15)  \
    ABS_TEST(a, b, tol);    \
  else                      \
    EXPECT_LE(std::abs(((a) - (b)) / (b)), tol);
