#ifndef LBM_LBM_BND_H
#define LBM_LBM_BND_H

#include <sfcmm_common.h>
#include "common/surface.h"


class LBMBndInterface {
 public:
  LBMBndInterface()          = default;
  virtual ~LBMBndInterface() = default;

  // deleted constructors not needed
  LBMBndInterface(const LBMBndInterface&) = delete;
  LBMBndInterface(LBMBndInterface&&)      = delete;
  auto operator=(const LBMBndInterface&) -> LBMBndInterface& = delete;
  auto operator=(LBMBndInterface&&) -> LBMBndInterface& = delete;

  virtual void init()  = 0;
  virtual void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) = 0;

 private:
};

template <LBMethodType LBTYPE, GBool TANGENTIALVELO>
class LBMBnd_wallBB;

template <Debug_Level DEBUG_LEVEL>
class LBMBndManager {
 public:
  LBMBndManager(const GInt dim, const GInt ndist) : m_dim(dim), m_ndist(ndist) {}
  virtual ~LBMBndManager() = default;

  // deleted constructors not needed
  LBMBndManager(const LBMBndManager&) = delete;
  LBMBndManager(LBMBndManager&&)      = delete;
  auto operator=(const LBMBndManager&) -> LBMBndManager& = delete;
  auto operator=(LBMBndManager&&) -> LBMBndManager& = delete;

  template <GInt NDIM>
  void addBndry(const BndryType bnd, const Surface<NDIM>& surf) {
    if constexpr(NDIM == 1) {
      if(m_ndist == 3) {
        addBndry<NDIM, 3>(bnd, surf);
      } else {
        TERMM(-1, "Invalid no of distributions!");
      }
    }
    if constexpr(NDIM == 2) {
      switch(m_ndist) {
        case 5:
          addBndry<NDIM, 5>(bnd, surf);
          break;
        case 9:
          addBndry<NDIM, 9>(bnd, surf);
          break;
        default:
          TERMM(-1, "Invalid number of distributions");
      }
    }
  }

  template <GInt NDIM, GInt NDIST>
  void addBndry(const BndryType bnd, const Surface<NDIM>& surf) {
    switch(bnd) {
      case BndryType::Wall_BounceBack:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<getLBMethodType<NDIM, NDIST>(), false>>(surf));
        break;
      case BndryType::Wall_BounceBack_TangentialVelocity:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<getLBMethodType<NDIM, NDIST>(), true>>(surf));
        break;
      default:
        TERMM(-1, "Invalid bndry Type!");
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold){
    for(const auto& bndry: m_bndrys){
      bndry->apply(f, fold);
    }
  }

 private:
  std::vector<std::unique_ptr<LBMBndInterface>> m_bndrys;

  GInt m_dim   = -1;
  GInt m_ndist = -1;
};

template <LBMethodType LBTYPE>
class LBMBndCell {
 public:
  LBMBndCell() = delete;
  //  LBMBndCell() = default;
  LBMBndCell(const GInt mappedCellId, const VectorD<dim(LBTYPE)>& normal) : m_mappedCellId(mappedCellId), m_normal(normal) {}
  virtual ~LBMBndCell() = default;

  virtual void init() {
    // calculate local bndry normal
  }

  [[nodiscard]] auto mapped() const -> GInt { return m_mappedCellId; }
  [[nodiscard]] auto normal() const -> const VectorD<dim(LBTYPE)>& { return m_normal; }

  // todo: fix me
  //  deleted constructors not needed
  //  LBMBndCell(const LBMBndCell&) = delete;
  //  LBMBndCell(LBMBndCell&&)      = delete;
  //  virtual auto operator=(const LBMBndCell&) -> LBMBndCell& = delete;
  //  virtual auto operator=(LBMBndCell&&) -> LBMBndCell& = delete;

 private:
  VectorD<dim(LBTYPE)> m_normal;
  GInt                 m_mappedCellId = -1;
};

template <LBMethodType LBTYPE, GBool TANGENTIALVELO>
class LBMBndCell_wallBB;

template <LBMethodType LBTYPE, GBool TANGENTIALVELO>
class LBMBnd_wallBB : public LBMBndInterface {
 public:
  LBMBnd_wallBB(const Surface<dim(LBTYPE)>& surf) {
    GInt surfId = 0;
    // todo: assign the surfaces to each bndry cell
    for(const GInt cellId : surf.getCellList()) {
      m_bndCells.emplace_back(cellId, surf.normal(surfId));
      ++surfId;
    }
    LBMBnd_wallBB<LBTYPE, TANGENTIALVELO>::init();
  }
  ~LBMBnd_wallBB() override = default;

  // deleted constructors not needed
  LBMBnd_wallBB(const LBMBnd_wallBB&) = delete;
  LBMBnd_wallBB(LBMBnd_wallBB&&)      = delete;
  auto operator=(const LBMBnd_wallBB&) -> LBMBnd_wallBB& = delete;
  auto operator=(LBMBnd_wallBB&&) -> LBMBnd_wallBB& = delete;

  void init() override {
    for(auto& bndCell : m_bndCells) {
      bndCell.init();
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) override {
    for(auto& bndCell : m_bndCells) {
      bndCell.apply(f, fold);
    }
  }

 private:
  VectorD<dim(LBTYPE)>                                   m_velocity;
  std::vector<LBMBndCell_wallBB<LBTYPE, TANGENTIALVELO>> m_bndCells;
};

template <LBMethodType LBTYPE, GBool TANGENTIALVELO>
class LBMBndCell_wallBB : public LBMBndCell<LBTYPE> {
  // class LBMBndCell_wallBB {
 public:
  LBMBndCell_wallBB() = delete;
  //  LBMBndCell_wallBB() = default;
  LBMBndCell_wallBB(const GInt cellId, const VectorD<dim(LBTYPE)>& _normal) : LBMBndCell<LBTYPE>(cellId, _normal) {}
  //  LBMBndCell_wallBB(const GInt cellId) {}
  virtual ~LBMBndCell_wallBB() = default;

  void init() {
    LBMBndCell<LBTYPE>::init();

    // precalculate weight and bndIndex
    for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
      if(inDirection<dim(LBTYPE)>(normal(), LBMethod<LBTYPE>::m_dirs[dist])) {
        m_bndIndex[m_noSetDists] = dist;
        ++m_noSetDists;
      }
    }
    //    m_weight
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) {
    for(GInt id = 0; id < m_noSetDists; ++id) {
      const GInt     dist          = m_bndIndex[id];
      GInt oppositeDist  = LBMethod<LBTYPE>::oppositeDist(dist);
      fold(mapped(), oppositeDist) = f(mapped(), dist);
    }
  }

  // todo: fix me
  //  deleted constructors not needed
  //  LBMBndCell_wallBB(const LBMBndCell_wallBB&) = delete;
  //  LBMBndCell_wallBB(LBMBndCell_wallBB&&)      = delete;
  //  auto operator=(const LBMBndCell_wallBB&) -> LBMBndCell_wallBB& = delete;
  //  auto operator=(LBMBndCell_wallBB&&) -> LBMBndCell_wallBB& = delete;

 private:
  using LBMBndCell<LBTYPE>::mapped;
  using LBMBndCell<LBTYPE>::normal;

  std::array<GDouble, noDists(LBTYPE)> m_weight;
  std::array<GInt, noDists(LBTYPE)>    m_bndIndex;
  GInt                                 m_noSetDists = 0;
};

#endif // LBM_LBM_BND_H
