#ifndef GRIDGENERATOR_FUNCTIONS_H
#define GRIDGENERATOR_FUNCTIONS_H
#include <bitset>
#include <gcem.hpp>
#include <iostream>
#include <json.h>
#include <map>
#include <sstream>
#include <unordered_set>
#include <sfcmm_common.h>
#include "common/term.h"

template <GInt NDIM, class T, class U>
static constexpr inline void assign(T& lhs, U& rhs) {
  for(int i = 0; i < NDIM; i++) {
    lhs[i] = static_cast<T>(rhs[i]);
  }
}

template <GInt NDIM, class T, class U>
static constexpr inline void fill(T& lhs, U value) {
  for(int i = 0; i < NDIM; i++) {
    lhs[i] = static_cast<T>(value);
  }
}

inline auto checkDuplicateIds(const std::vector<GInt>& ids) -> std::vector<GInt> {
  std::unordered_map<GInt, GInt> countMap;

  // Iterate over the vector and store the frequency of each element in map
  for(const auto& elem : ids) {
    const auto [it, success] = countMap.insert(std::pair<GInt, GInt>(elem, 1));
    if(!success) {
      // todo: using *it* produces a warning... (7/2021)
      countMap[elem]++;
    }
  }
  // Output for elements with more than count 1
  std::vector<GInt> duplicated;
  for(auto& it : countMap) {
    if(it.second > 1) {
      duplicated.emplace_back(it.first);
    }
  }
  return duplicated;
}

inline void removeDuplicates(std::vector<GInt>& id) {
  std::unordered_set<GInt> s;
  for(const auto i : id) {
    s.insert(i);
  }
  id.assign(s.begin(), s.end());
}

template <class T>
inline auto hasNAN(T& arrayToCheck) -> GInt {
  for(GInt id = 0; id < arrayToCheck.size(); ++id) {
    if(std::isnan(arrayToCheck[id])) {
      return id;
    }
  }
  return -1;
}
#endif // GRIDGENERATOR_FUNCTIONS_H
