#ifndef GRIDGENERATOR_BOUNDINGBOX_H
#define GRIDGENERATOR_BOUNDINGBOX_H

class BoundingBoxInterface {
 public:
  [[nodiscard]] virtual auto min(const GInt dir) const -> GDouble = 0;
  virtual auto               min(const GInt dir) -> GDouble&      = 0;

  [[nodiscard]] virtual inline auto max(const GInt dir) const -> GDouble = 0;
  virtual inline auto               max(const GInt dir) -> GDouble&      = 0;

  [[nodiscard]] virtual inline auto size() const -> GInt = 0;

  [[nodiscard]] virtual inline auto str() const -> GString = 0;
};


template <typename BoundingBoxType>
class BoundingBox : public BoundingBoxInterface {
 public:
  [[nodiscard]] inline auto min(const GInt dir) const -> GDouble override { return (*static_cast<const BoundingBoxType*>(this)).min[dir]; }

  inline auto min(const GInt dir) -> GDouble& override { return (*static_cast<BoundingBoxType*>(this)).min[dir]; }

  [[nodiscard]] inline auto max(const GInt dir) const -> GDouble override { return (*static_cast<const BoundingBoxType*>(this)).max[dir]; }

  inline auto max(const GInt dir) -> GDouble& override { return (*static_cast<BoundingBoxType*>(this)).max[dir]; }

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

template <GInt NDIM>
class BoundingBoxCT : public BoundingBox<BoundingBoxCT<NDIM>> {
 public:
  BoundingBoxCT() = default;
  BoundingBoxCT(const BoundingBoxDynamic& rhs);
  BoundingBoxCT(const BoundingBoxInterface& rhs);

  std::array<GDouble, NDIM> min{NAN_LIST<NDIM>()};
  std::array<GDouble, NDIM> max{NAN_LIST<NDIM>()};
  [[nodiscard]] inline auto size() const -> GInt override { return NDIM; }
};

class BoundingBoxReference : public BoundingBox<BoundingBoxReference> {
 public:
  GDouble*                  min = nullptr;
  GDouble*                  max = nullptr;
  [[nodiscard]] inline auto size() const -> GInt override { return ndim; }

  GInt ndim = 0;
};

class BoundingBoxDynamic : public BoundingBox<BoundingBoxDynamic> {
 public:
  BoundingBoxDynamic() = default;

  template <GInt NDIM>
  BoundingBoxDynamic(const BoundingBoxCT<NDIM>& rhs);

  void init(const GInt dim) {
    min.resize(dim);
    max.resize(dim);
    std::fill_n(min.begin(), dim, NAN);
    std::fill_n(max.begin(), dim, NAN);
  }

  [[nodiscard]] inline auto size() const -> GInt override { return min.size(); }


  std::vector<GDouble> min;
  std::vector<GDouble> max;
};


template <GInt NDIM>
BoundingBoxCT<NDIM>::BoundingBoxCT(const BoundingBoxDynamic& rhs) {
  for(GInt dir = 0; dir < NDIM; ++dir) {
    min[dir] = rhs.min[dir];
    max[dir] = rhs.max[dir];
  }
}

template <GInt NDIM>
BoundingBoxCT<NDIM>::BoundingBoxCT(const BoundingBoxInterface& rhs) {
  for(GInt dir = 0; dir < NDIM; ++dir) {
    min[dir] = rhs.min(dir);
    max[dir] = rhs.max(dir);
  }
}

template <GInt NDIM>
BoundingBoxDynamic::BoundingBoxDynamic(const BoundingBoxCT<NDIM>& rhs) {
  init(NDIM);
  for(GInt dir = 0; dir < NDIM; ++dir) {
    min[dir] = rhs.min[dir];
    max[dir] = rhs.max[dir];
  }
}

// template <GInt NDIM>
// inline auto VectorBBox(const BoundingBox<NDIM>& bbox) -> std::vector<GDouble> {
//   std::vector<GDouble> tmp;
//   for(GInt dir = 0; dir < NDIM; ++dir) {
//     tmp.emplace_back(bbox.min[dir]);
//     tmp.emplace_back(bbox.max[dir]);
//   }
//   return tmp;
// }


#endif // GRIDGENERATOR_BOUNDINGBOX_H
