// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_STRING_HELPER_H
#define SFCMM_STRING_HELPER_H
#include <bitset>
#include <gcem.hpp>
#include <iostream>
#include <map>
#include <sstream>
#include "common/sfcmm_types.h"

/// Get a stringstream from a given std::vector. Type needs to overload "<<".
/// \tparam LENGTH Length of the vector to be streamed.
/// \tparam T Type of the vector
/// \param in Vector to be stringstreamed
/// \return std::stringstream of the in vector.
template <class T>
static inline auto strStreamify(const std::vector<T>& in) -> std::stringstream {
  std::stringstream str;
  str << in[0];
  for(GUint i = 1; i < in.size(); i++) {
    str << " " << in[i];
  }
  return str;
}

/// Get a stringstream from a given std::vector. Type needs to overload "<<".
/// \tparam LENGTH Length of the vector to be streamed.
/// \tparam T Type of the vector
/// \param in Vector to be stringstreamed
/// \return std::stringstream of the in vector.
template <GInt LENGTH, class T>
static inline auto strStreamify(const std::vector<T>& in) -> std::stringstream {
  std::stringstream str;
  str << in[0];
  for(GInt i = 1; i < LENGTH; i++) {
    str << " " << in[i];
  }
  return str;
}

/// Get a stringstream from a given std::array. Type needs to overload "<<".
/// \tparam LENGTH Length of the vector to be streamed.
/// \tparam T Type of the array
/// \param in array to be stringstreamed
/// \return std::stringstream of the in array.
template <GInt LENGTH, class T>
static inline auto strStreamify(const std::array<T, LENGTH>& in) -> std::stringstream {
  std::stringstream str;
  str << in[0];
  for(GInt i = 1; i < LENGTH; i++) {
    str << " " << in[i];
  }
  return str;
}

/// Get a stringstream from a given std::array. Type needs to overload "<<".
/// \tparam LENGTH Length of the vector to be streamed.
/// \tparam T Type of the array
/// \param in array to be stringstreamed
/// \return std::stringstream of the in array.
template <GInt LENGTH>
static inline auto strStreamify(const VectorD<LENGTH>& in) -> std::stringstream {
  std::stringstream str;
  str << in[0];
  for(GInt i = 1; i < LENGTH; i++) {
    str << " " << in[i];
  }
  return str;
}

/// Convert an input vector to a string vector of the same size or as per the given optional argument size.
/// \tparam T Type of the vector.
/// \param in Input vector to be stringified.
/// \param size (Default=same as input vector) Provide the size if you want partial stringification.
/// \return Vector of strings of the input vector.
template <typename T>
static inline auto toStringVector(const std::vector<T>& in, GInt size = -1) -> std::vector<GString> {
  std::vector<GString> string_vector;

  if(size == -1) {
    size = in.size();
  }

  std::transform(in.begin(), in.begin() + size, std::back_inserter(string_vector), [](T b) -> GString { return std::to_string(b); });
  return string_vector;
}

/// Convert an input byte vector to a string vector of the same size or as per the given optional argument size. <std::byte version>
/// \param in Input byte vector to be stringified.
/// \param size (Default=same as input vector) Provide the size if you want partial stringification.
/// \return Vector of strings of the input byte vector.
template <>
inline auto toStringVector<std::byte>(const std::vector<std::byte>& in, GInt size) -> std::vector<GString> {
  std::vector<GString> string_vector;

  if(size == -1) {
    size = static_cast<GInt>(in.size());
  }

  std::transform(in.begin(), in.begin() + size, std::back_inserter(string_vector),
                 [](std::byte b) -> GString { return std::to_string(std::to_integer<GInt>(b)); });
  return string_vector;
}

/// Trim a string from the left.
/// \param s Input string
/// \return left trimmed string
static inline auto ltrim(std::string_view s) -> std::string_view {
  s.remove_prefix(std::distance(s.cbegin(), std::find_if(s.cbegin(), s.cend(), [](int c) { return std::isspace(c) == 0; })));

  return s;
}

/// Trim a string from the right.
/// \param s Input string
/// \return right trimmed string
static inline auto rtrim(std::string_view s) -> std::string_view {
  s.remove_suffix(std::distance(s.crbegin(), std::find_if(s.crbegin(), s.crend(), [](int c) { return std::isspace(c) == 0; })));

  return s;
}

/// Trim a string from both left and the right.
/// \param s Input string
/// \return Trimmed string
static inline auto trim(std::string_view s) -> std::string_view { return ltrim(rtrim(s)); }

/// Split a string by the given delimiter.
/// \tparam ContainerT Container type for the tokens produced by splitting the string by the delimiter
/// \param str String to be split.
/// \param tokens Container containing the string split into tokens.
/// \param delimiters Delimiter used for the splitting operation.
/// \param trimEmpty Should empty tokens be trimmed?
template <class ContainerT>
inline void tokenize(const std::string& str, ContainerT& tokens, const std::string& delimiters = " ", GBool trimEmpty = false) {
  std::string::size_type pos     = 0;
  std::string::size_type lastPos = 0;
  std::string::size_type length  = str.length();

  using value_type = typename ContainerT::value_type;
  using size_type  = typename ContainerT::size_type;

  while(lastPos < length + 1) {
    pos = str.find_first_of(delimiters, lastPos);
    if(pos == std::string::npos) {
      pos = length;
    }

    if(pos != lastPos || !trimEmpty) {
      tokens.push_back(value_type(str.data() + lastPos, static_cast<size_type>(pos) - lastPos));
    }

    lastPos = pos + 1;
  }
}

#endif // SFCMM_STRING_HELPER_H
