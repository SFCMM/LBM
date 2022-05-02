#ifndef COMMON_BASE64_H
#define COMMON_BASE64_H

#include "binary.h"
namespace base64 {
static constexpr GChar maxValue = 64;
// initialization value (0)
static constexpr GChar base64_zero{'A'};
// transposition table
static constexpr std::array<unsigned char, 65> encodeTable{"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"};

// Masks for accessing single bits
// 000001
static constexpr GInt mask_firstBit = 0x01;
// 000011
static constexpr GInt mask_twoBit = 0x03;
// 000111
static constexpr GInt mask_threeBit = 0x07;
// 001111
static constexpr GInt mask_fourBit = 0x0F;
// 011111
static constexpr GInt mask_fiveBit = 0x1F;
// 111111
static constexpr GInt mask_sixBit = 0x3F;
// 100000
static constexpr GInt mask_firstBitLE = 0x20;
// 110000
static constexpr GInt mask_twoBitLE = 0x30;
// 111000
static constexpr GInt mask_threeBitLE = 0x38;
// 111100
static constexpr GInt mask_fourBitLE = 0x3C;
// 111110
static constexpr GInt mask_fiveBitLE = 0x3E;

// tables to iterate access to single bits
static constexpr std::array<GInt, 7> masks   = {0, mask_firstBit, mask_twoBit, mask_threeBit, mask_fourBit, mask_fiveBit, mask_sixBit};
static constexpr std::array<GInt, 7> masksLE = {0,          mask_firstBitLE, mask_twoBitLE, mask_threeBitLE, mask_fourBitLE, mask_fiveBitLE,
                                                mask_sixBit};

/// Encode a single char to Base64
/// \param c Char to be encoded
/// \return Encoded char value in Base64
constexpr inline static unsigned char encodeChar(const unsigned char c) { return encodeTable[c]; }


/// Encode a variable of type T to Base64 chars. (Big Endian)
/// \tparam T Type of the variable to be encoded
/// \param c Value of the variable to be encoded
/// \return String of char values of the encoded type T variable
template <typename T>
inline static auto encode(const T c) -> GString {
  static constexpr GInt num_chars  = gcem::ceil(sizeof(T) * 8 / 6.0);
  static constexpr GInt shift      = (num_chars - 1) * 6;
  static constexpr GInt init_shift = sizeof(T) * 8 - shift;

  std::array<GUchar, num_chars> tmp{};
  T                             tmp_val = c;
  uint64_t                      tmp_int = *(static_cast<uint64_t*>(static_cast<void*>(&tmp_val)));
  tmp[0]                                = encodeTable[(tmp_int >> shift) & masks[init_shift]];
  for(int i = 1; i < num_chars; ++i) {
    tmp[i] = encodeTable[(tmp_int >> (shift - i * 6)) & masks[6]];
  }
  return {std::begin(tmp), std::end(tmp)};
}

/// Encode an array of type T variables to Base64  (Big Endian)
/// \tparam T Type of the variables to be encoded.
/// \tparam length Length of the array to be encoded.
/// \param c Values to be encoded.
/// \return String of the encoded array.
template <typename T, GInt length>
inline static auto encode(T* c) -> GString {
  static constexpr GInt num_chars = gcem::ceil(sizeof(T) * 8 * length / 6.0);

  std::array<GUchar, num_chars> encoded_base64{};
  std::bitset<num_chars * 6>    mem{};

  auto* char_wise = static_cast<GUchar*>(static_cast<void*>(&c[0]));

  for(GInt i = 0; i < length; ++i) {
    for(GUint byte = 0; byte < sizeof(T); ++byte) {
      auto tmp_bitset = std::bitset<8>(char_wise[i * sizeof(T) + byte]);
      for(GInt bit = 0; bit < 8; ++bit) {
        mem[(length - i - 1) * sizeof(T) * 8 + byte * 8 + bit] = tmp_bitset[bit];
      }
    }
  }

  for(GInt i = 0; i < num_chars; ++i) {
    const GInt num = mem[i * 6] + mem[i * 6 + 1] * 2 + mem[i * 6 + 2] * 4 + mem[i * 6 + 3] * 8 + mem[i * 6 + 4] * 16 + mem[i * 6 + 5] * 32;
    encoded_base64[num_chars - i - 1] = encodeTable[num];
  }

  return {std::begin(encoded_base64), std::end(encoded_base64)};
}

/// Encode a variable of type T to Base64 chars. (Little Endian)
/// \tparam T Type of the variable to be encoded
/// \param c Value of the variable to be encoded
/// \return String of char values of the encoded type T variable
template <typename T>
inline static auto encodeLE(const T c) -> GString {
  static constexpr GInt num_chars = gcem::ceil(sizeof(T) * 8 / 6.0);
  static constexpr GInt shift     = sizeof(T) * 8 - 6;
  static constexpr GInt end_shift = num_chars * 6 - sizeof(T) * 8;

  std::array<GUchar, num_chars> tmp{};
  T                             tmp_val = binary::getSwappedEndian(c);
  uint64_t                      tmp_int = *(static_cast<uint64_t*>(static_cast<void*>(&tmp_val)));
  for(int i = 0; i < num_chars - 1; ++i) {
    tmp[i] = encodeTable[(tmp_int >> (shift - i * 6)) & masksLE[6]];
  }
  tmp[num_chars - 1] = encodeTable[(tmp_int << end_shift) & masksLE[6 - end_shift]];
  return {std::begin(tmp), std::end(tmp)};
}

/// Encode an array of type T variables to Base64  (Little Endian)
/// \tparam T Type of the variables to be encoded.
/// \tparam length Length of the array to be encoded.
/// \tparam shifted Pad at the beginning to align with byte boundary.
/// \param c Values to be encoded.
/// \return String of the encoded array.
// template <typename T, GInt length, GInt shifted = 0>
// inline static auto encodeLE(const T* c) -> GString {
//   static constexpr GInt num_chars = gcem::ceil((sizeof(T) * 8 * length - shifted) / 6.0);
//
//   std::array<T, length>         swapped_endian{};
//   std::array<GUchar, num_chars> encoded_base64{};
//   std::bitset<num_chars * 6>    mem{};
//   for(GInt i = 0; i < length; ++i) {
//     swapped_endian[i] = binary::getSwappedEndian(c[i]);
//   }
//
//   auto* char_wise = static_cast<GUchar*>(static_cast<void*>(&swapped_endian[0]));
//
//   for(GInt i = 0; i < length; ++i) {
//     for(GUint byte = 0; byte < sizeof(T); ++byte) {
//       auto tmp_bitset = std::bitset<8>(char_wise[i * sizeof(T) + byte]);
//       for(GInt bit = 0; bit < 8; ++bit) {
//         const GInt index = num_chars * 6 - (i + 1) * sizeof(T) * 8 + byte * 8 + bit + shifted;
//         if(shifted == 0 || index < num_chars * 6) {
//           mem[index] = tmp_bitset[bit];
//         } else if(shifted > 0 && tmp_bitset[bit]) {
//           std::cerr << "WARNING: could be wrong" << std::endl; // todo: remove
//           TERMM(-1, "writing is incorrect");
//         }
//       }
//     }
//   }
//
//   for(GInt i = 0; i < num_chars; ++i) {
//     const GInt num = mem[i * 6] + mem[i * 6 + 1] * 2 + mem[i * 6 + 2] * 4 + mem[i * 6 + 3] * 8 + mem[i * 6 + 4] * 16 + mem[i * 6 + 5] *
//     32; encoded_base64[num_chars - i - 1] = encodeTable[num];
//   }
//
//   return {std::begin(encoded_base64), std::end(encoded_base64)};
// }

// todo: unify with the compile-time constant version!

/// Encode an array of type T variables to Base64  (Little Endian) (non compile-time constant version)
/// \tparam T Type of the variables to be encoded.
/// \param length Length of the array to be encoded.
/// \param c Values to be encoded.
/// \return String of the encoded array.
template <typename T>
inline static auto encodeLE(const T* c, const GInt length) -> GString {
  // Nothing to do for length 0
  if(length == 0) {
    std::cerr << "ERROR: Invalid call to encodeLE() with length = 0" << std::endl;
    std::exit(-1);
    // return "";
  }

  const auto num_chars = static_cast<GUint>(gcem::ceil((sizeof(T) * 8 * length) / 6.0));

  std::vector<T>      swapped_endian(length);
  std::vector<GUchar> encoded_base64(num_chars);
  std::vector<GBool>  mem(num_chars * 6);

  // perform the actual endian swap
  for(GInt i = 0; i < length; ++i) {
    swapped_endian[i] = binary::getSwappedEndian(c[i]);
  }

  // perform byte alignment and switch to bit representation
  auto* char_wise = static_cast<GUchar*>(static_cast<void*>(swapped_endian.data()));
  for(GInt i = 0; i < length; ++i) {
    for(GUint byte = 0; byte < sizeof(T); ++byte) {
      auto tmp_bitset = std::bitset<8>(char_wise[i * sizeof(T) + byte]);
      for(GInt bit = 0; bit < 8; ++bit) {
        const GInt index = num_chars * 6 - (i + 1) * sizeof(T) * 8 + byte * 8 + bit;
        mem[index]       = tmp_bitset[bit];
      }
    }
  }


  for(GInt i = 0; i < num_chars; ++i) {
    // convert 6-bit representation to a number
    const GInt num = mem[i * 6] + mem[i * 6 + 1] * 2 + mem[i * 6 + 2] * 4 + mem[i * 6 + 3] * 8 + mem[i * 6 + 4] * 16 + mem[i * 6 + 5] * 32;
    // reverse mem and encode to char
    encoded_base64[num_chars - i - 1] = encodeTable[num];
  }

  return {std::begin(encoded_base64), std::end(encoded_base64)};
}

/// Encode an array of type T variables to Base64  (Little Endian) (non compile-time constant version)
/// \tparam T Type of the variables to be encoded.
/// \tparam shifted Number of bits to pad at the beginning to align with byte boundary.
/// \param length Length of the array to be encoded.
/// \param c Values to be encoded.
/// \return String of the encoded array.
template <typename T, typename U>
inline static auto encodeLE_header(const T* c, const U length) -> GString {
  // Nothing to do for length 0
  if(length == 0) {
    std::cerr << "ERROR: Invalid call to encodeLE() with length = 0" << std::endl;
    std::exit(-1);
    // return "";
  }

  const auto  num_chars      = static_cast<GUint>(gcem::ceil((sizeof(T) * 8 * length + sizeof(U) * 8) / 6.0));
  const GUint actual_bit     = sizeof(T) * 8 * length + sizeof(U) * 8;
  const GUint allocated_bits = num_chars * 6;

  U                   swapped_length = 0;
  std::vector<T>      swapped_endian(length);
  std::vector<GUchar> encoded_base64(num_chars);
  std::vector<GBool>  mem(num_chars * 6);

  // perform the actual endian swap
  swapped_length = binary::getSwappedEndian(length * binary::BYTE_SIZE);
  for(GInt i = 0; i < length; ++i) {
    swapped_endian[i] = binary::getSwappedEndian(c[i]);
  }

  // perform byte alignment and switch to bit representation
  // todo: remove one loop
  GUint count     = allocated_bits - actual_bit;
  auto* char_wise = static_cast<GUchar*>(static_cast<void*>(swapped_endian.data()));
  for(GInt i = 0; i < length; ++i) {
    for(GUint byte = 0; byte < sizeof(T); ++byte) {
      // access bytewise the swapped memory buffer
      auto tmp_bitset = std::bitset<8>(char_wise[(length - i - 1) * sizeof(T) + byte]);
      for(GInt bit = 0; bit < 8; ++bit) {
        mem.at(count) = tmp_bitset[bit];
        ++count;
      }
    }
  }

  // add length information to the beginning
  // todo: remove loop
  auto tmp_bitset = std::bitset<sizeof(U) * 8>(swapped_length);
  for(GUint bit = 0; bit < sizeof(U) * 8; ++bit) {
    mem.at(count) = tmp_bitset[bit];
    ++count;
  }


  for(GUint i = 0; i < num_chars; ++i) {
    // convert 6-bit representation to a number
    const GInt num = mem[i * 6] + mem[i * 6 + 1] * 2 + mem[i * 6 + 2] * 4 + mem[i * 6 + 3] * 8 + mem[i * 6 + 4] * 16 + mem[i * 6 + 5] * 32;
    // reverse mem and encode to char
    encoded_base64[num_chars - i - 1] = encodeTable[num];
  }

  return {std::begin(encoded_base64), std::end(encoded_base64)};
}

} // namespace base64

#endif // COMMON_BASE64_H
