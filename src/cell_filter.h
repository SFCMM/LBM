#ifndef LBM_CELL_FILTER_H
#define LBM_CELL_FILTER_H

#include <json.h>
#include <sfcmm_common.h>
#include "cartesiangrid_base.h"

class FilterBase {
 public:
  virtual ~FilterBase() = default;
  [[nodiscard]] virtual auto eval(const GInt /*id*/) const -> GBool { return true; };
};

template <GInt NDIM>
class CellFilterBase : public FilterBase {
 public:
  CellFilterBase(CartesianGridData<NDIM> grid) : m_grid(grid) {}

  [[nodiscard]] auto eval(const GInt /*cellId*/) const -> GBool override { return true; };

  auto grid() const -> const CartesianGridData<NDIM>& { return m_grid; }

 private:
  CartesianGridData<NDIM> m_grid;
};

template <GInt NDIM>
class LeafCellFilter : public CellFilterBase<NDIM> {
 public:
  LeafCellFilter(CartesianGridData<NDIM> grid) : CellFilterBase<NDIM>(grid) {}

  [[nodiscard]] auto eval(const GInt cellId) const -> GBool override { return CellFilterBase<NDIM>::grid().isLeaf(cellId); }
};

template <GInt NDIM>
class LevelCellFilter : public CellFilterBase<NDIM> {
 public:
  LevelCellFilter(CartesianGridData<NDIM> grid, const GInt lvl) : CellFilterBase<NDIM>(grid), m_targetLvl(lvl) {}

  [[nodiscard]] auto eval(const GInt cellId) const -> GBool override {
    return std::to_integer<GInt>(CellFilterBase<NDIM>::grid().level(cellId)) == m_targetLvl;
  }

 private:
  GInt m_targetLvl;
};

template <GInt NDIM>
class AllFilter : public FilterBase {
 public:
  AllFilter() = default;

  [[nodiscard]] auto eval(const GInt /*cellId*/) const -> GBool override { return true; }
};

template <GInt NDIM>
class CellFilterManager {
 private:
  using json = nlohmann::json;

 public:
  CellFilterManager() { m_filter.emplace_back(std::make_unique<AllFilter<NDIM>>()); }

  CellFilterManager(const json& /*filterConf*/) { m_filter.emplace_back(std::make_unique<AllFilter<NDIM>>()); }

  CellFilterManager(const json& filterConf, CartesianGridData<NDIM> grid) {
    if(config::has_config_value(filterConf, "cellFilter")) {
      json cellFilterConf = filterConf["cellFilter"];
      if(cellFilterConf.is_array()) {
        for(const auto& filter : std::vector<GString>(filterConf["cellFilter"])) {
          addFilter(filter, filterConf, grid);
        }
        // todo: add test
        TERMM(-1, "untested");
      } else {
        addFilter(filterConf["cellFilter"], filterConf, grid);
      }
    } else {
      TERMM(-1, "Invalid output configuration missing \"cellFilter\" key");
    }
  }

  auto addFilter(const GString& filter, const json& filterConf, CartesianGridData<NDIM> grid) {
    // todo: replace with switch
    if(filter == "highestLvl") {
      m_filter.emplace_back(std::make_unique<LevelCellFilter<NDIM>>(grid, grid.currentHighestLvl()));
    } else if(filter == "lowestLvl" || filter == "partitionLvl") {
      m_filter.emplace_back(std::make_unique<LevelCellFilter<NDIM>>(grid, grid.partitionLvl()));
    } else if(filter == "leafCells") {
      m_filter.emplace_back(std::make_unique<LeafCellFilter<NDIM>>(grid));
    } else if(filter == "targetLvl") {
      const GInt outputLvl = config::required_config_value<GInt>(filterConf, "outputLvl");
      if(outputLvl < grid.partitionLvl()) {
        TERMM(-1, "Outputting a lvl below the partition lvl is not possible!");
      }
      m_filter.emplace_back(std::make_unique<LevelCellFilter<NDIM>>(grid, outputLvl));
    } else {
      TERMM(-1, "Unknown output filter " + static_cast<GString>(filter));
    }
  }

  [[nodiscard]] auto eval(const GInt cellId) const -> GBool {
    for(const auto& filter : m_filter) {
      if(!filter->eval(cellId)) {
        return false;
      }
    }
    return true;
  }

 private:
  std::vector<std::unique_ptr<FilterBase>> m_filter;
};

#endif // LBM_CELL_FILTER_H
