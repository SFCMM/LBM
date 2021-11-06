// SPDX-License-Identifier: BSD-3-Clause

#ifndef COMMON_BINARY_H
#define COMMON_BINARY_H

#include <bitset>

namespace binary {
// number of bits in a byte
static constexpr GShort BYTE_SIZE = 8;

/// Convert a value to the binary representation
/// \tparam T Type of the value to be converted.
/// \param num Value to be converted.
/// \return Bit representation of the value.
template <typename T>
static constexpr inline auto convert(const T num) -> std::bitset<sizeof(T) * BYTE_SIZE> {
  T tmp = num;
  return std::bitset<sizeof(T) * BYTE_SIZE>(*(static_cast<uint64_t*>(static_cast<void*>(&tmp))));
}
#if defined(GCC_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstack-protector"
#endif
/// Swap the endianess of a variable.
/// \tparam T Type to endianess swapped.
/// \param val Value to endianess swapped.
/// \return Return value but with swapped endianess.
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

/// Swap the endianess of a variable. (constant variable version)
/// \tparam T Type to endianess swapped.
/// \param val Value to endianess swapped.
/// \return Return value but with swapped endianess.
template <typename T>
static constexpr inline auto getSwappedEndian(const T val) {
  T tmp = val;
  swapEndian(tmp);
  return tmp;
}

/// Special version of endian swapping for 16bit variables.
/// \param value 16bit variable.
/// \return Swapped 16bit variable (in-place)
template <>
constexpr inline void swapEndian<std::uint16_t>(std::uint16_t& value) {
  value = (value >> BYTE_SIZE) | (value << BYTE_SIZE);
}

/// Special version of endian swapping for 32bit variables.
/// \param value 32bit variable.
/// \return Swapped 32bit variable (in-place)
template <>
constexpr inline void swapEndian<std::uint32_t>(std::uint32_t& value) {
  std::uint32_t tmp = ((value << BYTE_SIZE) & 0xFF00FF00) | ((value >> BYTE_SIZE) & 0xFF00FF);
  value             = (tmp << 2 * BYTE_SIZE) | (tmp >> 2 * BYTE_SIZE);
}

/// Special version of endian swapping for 64bit variables.
/// \param value 64bit variable.
/// \return Swapped 64bit variable (in-place)
template <>
constexpr inline void swapEndian<std::uint64_t>(std::uint64_t& value) {
  value = ((value & 0x00000000FFFFFFFFull) << 4 * BYTE_SIZE) | ((value & 0xFFFFFFFF00000000ull) >> 4 * BYTE_SIZE);
  value = ((value & 0x0000FFFF0000FFFFull) << 2 * BYTE_SIZE) | ((value & 0xFFFF0000FFFF0000ull) >> 2 * BYTE_SIZE);
  value = ((value & 0x00FF00FF00FF00FFull) << BYTE_SIZE) | ((value & 0xFF00FF00FF00FF00ull) >> BYTE_SIZE);
}

/// Convert a variable into endianess swapped binary representation.
/// \tparam T Type of the variable to endianess swapped.
/// \param num Value to be swapped.
/// \return Bit representation of the swapped variable.
template <typename T>
static constexpr inline auto convertSwap(const T num) -> std::bitset<sizeof(T) * BYTE_SIZE> {
  T tmp = num;
  swapEndian(tmp);
  return convert(tmp);
}

} // namespace binary

#endif // COMMON_BINARY_H
