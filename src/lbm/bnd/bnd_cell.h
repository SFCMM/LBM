#ifndef LBM_BND_CELL_H
#define LBM_BND_CELL_H

#include "lbm/constants.h"
#include "sfcmm_common.h"

template <LBMethodType LBTYPE>
class LBMBndCell {
 public:
  LBMBndCell(const GInt mappedCellId, const VectorD<dim(LBTYPE)>& normal) : m_normal(normal), m_mappedCellId(mappedCellId) {
    m_outsideDir.fill(INVALID_DIR);
  }
  virtual ~LBMBndCell() = default;

  LBMBndCell(const LBMBndCell&)                            = default;
  LBMBndCell(LBMBndCell&&) noexcept                        = default;
  virtual auto operator=(const LBMBndCell&) -> LBMBndCell& = delete;
  virtual auto operator=(LBMBndCell&&) -> LBMBndCell&      = delete;

  virtual void init() {
    // precalculate weight and bndIndex
    for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
      // determine if the dist points into the outside normal direction of bndry
      if(inDirection<dim(LBTYPE)>(normal(), LBMethod<LBTYPE>::m_dirs[dist])) {
        m_outsideDir[dist] = dist;
      }
    }
  }

  /** Map the bndCell to the cellId of the local cell storage
   * @return The mapped cellId of the local cell storage for this bndCell
   */
  [[nodiscard]] auto mapped() const -> GInt { return m_mappedCellId; }

  /** Return the normal surface vector pointing to the outside of the geometry
   * @return Normal vector pointing outside
   */
  [[nodiscard]] auto normal() const -> const VectorD<dim(LBTYPE)>& { return m_normal; }

  /** Return the orientation pointing to the outside of the geometry
   * @return Orientation pointing outside
   */
  [[nodiscard]] auto outsideDir(const GInt dist) const -> GInt { return m_outsideDir[dist]; }

 private:
  VectorD<dim(LBTYPE)>              m_normal;
  GInt                              m_mappedCellId = -1;
  std::array<GInt, noDists(LBTYPE)> m_outsideDir;
};

#endif // LBM_BND_CELL_H
