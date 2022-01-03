#ifndef LBM_RANDOM_H
#define LBM_RANDOM_H

#include <sfcmm_common.h>

class RandomGenerator {
 public:
  RandomGenerator() : m_generator(123456) { logger << "WARNING: No seed set" << std::endl; };

  template <GInt NDIM>
  inline auto randomPos(const VectorD<NDIM>& lowerBound, const VectorD<NDIM>& upperBound) -> VectorD<NDIM> {
    VectorD<NDIM> position;
    for(GInt dir = 0; dir < NDIM; ++dir) {
      position[dir] = m_generator.double_value(lowerBound[dir], upperBound[dir]);
    }
    return position;
  }

 private:
  randxor m_generator;
};

#endif // LBM_RANDOM_H
