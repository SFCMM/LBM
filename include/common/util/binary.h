// SPDX-License-Identifier: BSD-3-Clause

#ifndef COMMON_BINARY_H
#define COMMON_BINARY_H

#include <bitset>

namespace binary {
static constexpr GShort BYTE_SIZE = 8;

template <typename T>
static constexpr inline auto convert(const T num) -> std::bitset<sizeof(T) * BYTE_SIZE> {
  T tmp = num;
  return std::bitset<sizeof(T) * BYTE_SIZE>(*(static_cast<uint64_t*>(static_cast<void*>(&tmp))));
}
#if defined(GCC_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstack-protector"
#endif
template <typename T>
static constexpr inline void swapEndian(T& val) {
  union U {
    T                                   val;
    std::array<std::uint8_t, sizeof(T)> raw;
  } src, dst;

  src.val = val;
  std::reverse_copy(src.raw.begin(), src.raw.end(), dst.raw.begin());
  val = dst.val;
}
#if defined(GCC_COMPILER)
#pragma GCC diagnostic pop
#endif

template <typename T>
static constexpr inline auto getSwappedEndian(const T val) {
  T tmp = val;
  swapEndian(tmp);
  return tmp;
}

template <>
constexpr inline void swapEndian<std::uint16_t>(std::uint16_t& value) {
  value = (value >> BYTE_SIZE) | (value << BYTE_SIZE);
}

template <>
constexpr inline void swapEndian<std::uint32_t>(std::uint32_t& value) {
  std::uint32_t tmp = ((value << BYTE_SIZE) & 0xFF00FF00) | ((value >> BYTE_SIZE) & 0xFF00FF);
  value             = (tmp << 2 * BYTE_SIZE) | (tmp >> 2 * BYTE_SIZE);
}

template <>
constexpr inline void swapEndian<std::uint64_t>(std::uint64_t& value) {
  value = ((value & 0x00000000FFFFFFFFull) << 4 * BYTE_SIZE) | ((value & 0xFFFFFFFF00000000ull) >> 4 * BYTE_SIZE);
  value = ((value & 0x0000FFFF0000FFFFull) << 2 * BYTE_SIZE) | ((value & 0xFFFF0000FFFF0000ull) >> 2 * BYTE_SIZE);
  value = ((value & 0x00FF00FF00FF00FFull) << BYTE_SIZE) | ((value & 0xFF00FF00FF00FF00ull) >> BYTE_SIZE);
}

template <typename T>
static constexpr inline auto convertSwap(const T num) -> std::bitset<sizeof(T) * BYTE_SIZE> {
  T tmp = num;
  swapEndian(tmp);
  return convert(tmp);
}

} // namespace binary

#endif // COMMON_BINARY_H
