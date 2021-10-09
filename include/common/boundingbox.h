#ifndef GRIDGENERATOR_BOUNDINGBOX_H
#define GRIDGENERATOR_BOUNDINGBOX_H

#include "common/constants.h"
#include "common/sfcmm_types.h"

/// Interface for accessing bounding boxes.
class BoundingBoxInterface {
 public:
  /// Return the minimum values of the bounding box (const)
  /// \param dir Direction of the value
  /// \return min[dir]
  [[nodiscard]] virtual auto min(const GInt dir) const -> GDouble = 0;

  /// Return the minimum values of the bounding box
  /// \param dir Direction of the value
  /// \return min[dir]
  virtual auto min(const GInt dir) -> GDouble& = 0;

  /// Return the maximum values of the bounding box (const)
  /// \param dir Direction of the value
  /// \return max[dir]
  [[nodiscard]] virtual inline auto max(const GInt dir) const -> GDouble = 0;

  /// Return the maximum values of the bounding box
  /// \param dir Direction of the value
  /// \return max[dir]
  virtual inline auto max(const GInt dir) -> GDouble& = 0;

  /// Number of values in the bounding box (most of the time size = NDIM)
  /// \return Number of directions.
  [[nodiscard]] virtual inline auto size() const -> GInt = 0;

  /// Stringify the bounding box.
  /// \return String version of all the values.
  [[nodiscard]] virtual inline auto str() const -> GString = 0;
};


/// CRTP Base-Class for the bounding box memory classes below.
/// \tparam BoundingBoxType CRTP for the container classes below.
template <typename BoundingBoxType>
class BoundingBox : public BoundingBoxInterface {
 public:
  [[nodiscard]] inline auto min(const GInt dir) const -> GDouble override {
    return (*static_cast<const BoundingBoxType*>(this)).m_min[dir];
  }

  inline auto min(const GInt dir) -> GDouble& override { return (*static_cast<BoundingBoxType*>(this)).m_min[dir]; }

  [[nodiscard]] inline auto max(const GInt dir) const -> GDouble override {
    return (*static_cast<const BoundingBoxType*>(this)).m_max[dir];
  }

  inline auto max(const GInt dir) -> GDouble& override { return (*static_cast<BoundingBoxType*>(this)).m_max[dir]; }

  [[nodiscard]] inline auto size() const -> GInt override { return static_cast<const BoundingBoxType*>(this)->size(); }

  [[nodiscard]] inline auto str() const -> GString override {
    GString text = "{";
    for(GInt i = 0; i < size(); ++i) {
      text += "[" + std::to_string(min(i)) + ", " + std::to_string(max(i)) + "]";
    }
    text += "}";
    return text;
  }
};

class BoundingBoxDynamic;

/// Compile time constant version of the bounding box class.
/// \tparam NDIM Number of directions of the bounding box.
template <GInt NDIM>
class BoundingBoxCT : public BoundingBox<BoundingBoxCT<NDIM>> {
  friend BoundingBox<BoundingBoxCT<NDIM>>;

 public:
  BoundingBoxCT() = default;
  /// Initialize from dynamic bounding box
  /// \param rhs Dynamic Bounding Box type to copy from
  BoundingBoxCT(const BoundingBoxDynamic& rhs);

  /// Initialize from any bounding box
  /// \param rhs Any bounding box type
  BoundingBoxCT(const BoundingBoxInterface& rhs);

  [[nodiscard]] inline auto size() const -> GInt override { return NDIM; }

 private:
  std::array<GDouble, NDIM> m_min{NAN_LIST<NDIM>()};
  std::array<GDouble, NDIM> m_max{NAN_LIST<NDIM>()};
};

/// Reference type for bounding boxes.
class BoundingBoxReference : public BoundingBox<BoundingBoxReference> {
  friend BoundingBox<BoundingBoxReference>;

 public:
  [[nodiscard]] inline auto size() const -> GInt override { return ndim; }

 private:
  GInt ndim = 0;

  GDouble* min = nullptr;
  GDouble* max = nullptr;
};

/// Dynamic memory version of the bounding box type.
class BoundingBoxDynamic : public BoundingBox<BoundingBoxDynamic> {
  friend BoundingBox<BoundingBoxDynamic>;

 public:
  // todo: add constructor with dimension
  BoundingBoxDynamic() = default;

  /// Initialize from a compile time constant bounding box
  /// \param rhs Compile time constant bounding box type
  template <GInt NDIM>
  BoundingBoxDynamic(const BoundingBoxCT<NDIM>& rhs);

  /// initialize dynamic type
  /// \param dim Dimension of the initialization.
  // todo: add init to interface
  void init(const GInt dim) {
    m_min.resize(dim);
    m_max.resize(dim);
    std::fill_n(m_min.begin(), dim, NAN);
    std::fill_n(m_max.begin(), dim, NAN);
  }

  [[nodiscard]] inline auto size() const -> GInt override { return static_cast<GInt>(m_min.size()); }

 private:
  std::vector<GDouble> m_min;
  std::vector<GDouble> m_max;
};


template <GInt NDIM>
BoundingBoxCT<NDIM>::BoundingBoxCT(const BoundingBoxDynamic& rhs) {
  for(GInt dir = 0; dir < NDIM; ++dir) {
    m_min[dir] = rhs.min(dir);
    m_max[dir] = rhs.max(dir);
  }
}

template <GInt NDIM>
BoundingBoxCT<NDIM>::BoundingBoxCT(const BoundingBoxInterface& rhs) {
  for(GInt dir = 0; dir < NDIM; ++dir) {
    m_min[dir] = rhs.min(dir);
    m_max[dir] = rhs.max(dir);
  }
}

template <GInt NDIM>
BoundingBoxDynamic::BoundingBoxDynamic(const BoundingBoxCT<NDIM>& rhs) {
  init(NDIM);
  for(GInt dir = 0; dir < NDIM; ++dir) {
    m_min[dir] = rhs.min(dir);
    m_max[dir] = rhs.max(dir);
  }
}


#endif // GRIDGENERATOR_BOUNDINGBOX_H
