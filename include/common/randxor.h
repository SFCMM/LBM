// SPDX-License-Identifier: BSD-3-Clause

#ifndef COMMON_RANDXOR_H
#define COMMON_RANDXOR_H

#include "common/sfcmm_types.h"

/// State of splitmix64
struct splitmix64_state {
  uint64_t s;
};

/// State of xoshiro256+
struct xoshiro256p_state {
  std::array<uint64_t, 4> s = {1, 2, 3, 4};
};

/// splitmix is used for initialization
/// \param state Splitmix state
/// \return Pseudorandom value of the splitmix64 function
inline uint64_t splitmix64(splitmix64_state* state) {
  uint64_t result = (state->s += 0x9E3779B97f4A7C15);
  result          = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
  result          = (result ^ (result >> 27)) * 0x94D049BB133111EB;
  return result ^ (result >> 31);
}

///
/// \param x
/// \param k
/// \return
inline uint64_t rol64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

/// Advance the state of the xoshiro prng
/// \param state State of the xoshiro prng
/// \return Pseudo Random value based on the given state function
inline uint64_t xoshiro256p(xoshiro256p_state* state) {
  uint64_t*      s      = state->s.data();
  const uint64_t result = s[0] + s[3];
  const uint64_t t      = s[1] << 17;

  s[2] ^= s[0];
  s[3] ^= s[1];
  s[1] ^= s[2];
  s[0] ^= s[3];

  s[2] ^= t;
  s[3] = rol64(s[3], 45);

  return result;
}

/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */
inline void jump(xoshiro256p_state* state) {
  static const uint64_t JUMP[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c};

  uint64_t* s  = state->s.data();
  uint64_t  s0 = 0;
  uint64_t  s1 = 0;
  uint64_t  s2 = 0;
  uint64_t  s3 = 0;
  for(uint64_t i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
    for(int b = 0; b < 64; b++) {
      if(JUMP[i] & UINT64_C(1) << b) {
        s0 ^= s[0];
        s1 ^= s[1];
        s2 ^= s[2];
        s3 ^= s[3];
      }
      xoshiro256p(state);
    }

  s[0] = s0;
  s[1] = s1;
  s[2] = s2;
  s[3] = s3;
}


/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^64 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations. */
inline void long_jump(xoshiro256p_state* state) {
  static const uint64_t LONG_JUMP[] = {0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635};

  uint64_t* s  = state->s.data();
  uint64_t  s0 = 0;
  uint64_t  s1 = 0;
  uint64_t  s2 = 0;
  uint64_t  s3 = 0;
  for(uint64_t i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
    for(int b = 0; b < 64; b++) {
      if(LONG_JUMP[i] & UINT64_C(1) << b) {
        s0 ^= s[0];
        s1 ^= s[1];
        s2 ^= s[2];
        s3 ^= s[3];
      }
      xoshiro256p(state);
    }

  s[0] = s0;
  s[1] = s1;
  s[2] = s2;
  s[3] = s3;
}

class randxor {
 public:
  explicit randxor(uint64_t seed) : initState{seed} {
    //    initState          = {seed};
    internalState.s[0] = splitmix64(&initState);
    internalState.s[1] = splitmix64(&initState);
    internalState.s[2] = splitmix64(&initState);
    internalState.s[3] = splitmix64(&initState);
  }

  /// Generate random integer value
  /// \return random 64bit integer value
  inline auto uinteger_value() -> uint64_t { return xoshiro256p(&internalState); }

  /// Generate random integer in [lowerbound, upperbound)
  /// \param lower_bound lower bound
  /// \param upper_bound upper bound
  /// \return Random unsigned integer in [lowerbound, upperbound)
  inline auto uinteger_value(uint32_t lower_bound, uint32_t upper_bound) -> uint64_t {
    assert(upper_bound > lower_bound);

    return ((static_cast<uint64_t>(static_cast<uint32_t>(uinteger_value())) * static_cast<uint64_t>(upper_bound - lower_bound)) >> 32)
           + static_cast<uint64_t>(lower_bound);
  }


  /// Generate random double value
  /// \return random value in [0, 1]
  inline auto double_value() -> GDouble {
    static constexpr uint64_t DOUBLE_MASK = UINT64_C(0x3FF) << 52;

    union {
      uint64_t int_value;
      GDouble  val;
    };

    int_value = (xoshiro256p(&internalState) >> 12) | DOUBLE_MASK;

    return val - 1.0;
  }

  /// Generate random double value in [lowerbound, upperbound)
  /// \param lower_bound lower bound
  /// \param upper_bound upper bound
  /// \return Random double value in [lower_bound, upper_bound)
  inline double double_value(double lower_bound, double upper_bound) {
    return (double_value() * (upper_bound - lower_bound)) + lower_bound;
  }

 private:
  splitmix64_state  initState;
  xoshiro256p_state internalState;
};


#endif // COMMON_RANDXOR_H
