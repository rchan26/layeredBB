#ifndef DOUBLE_COMPARISON
#define DOUBLE_COMPARISON

#include <cmath>
#include <limits>

bool almostEqual(const double &a,
                 const double &b,
                 const double &epsilon)
{
  return (std::fabs(a-b) <= epsilon);
}

bool approximatelyEqualAbsRel(const double &a,
                              const double &b,
                              const double &absEpsilon,
                              const double &relEpsilon)
{
  const double diff = std::fabs(a-b);
  if (diff <= absEpsilon) {
    // Check if the numbers are really close -- needed when comparing numbers near zero.
    return true;
  } else {
    // Otherwise fall back to Knuth's algorithm for relative equality
    return (diff <= (std::max(std::abs(a),std::abs(b)) * relEpsilon));
  }
}

#endif