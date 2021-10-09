#ifndef GRIDGENERATOR_GRIDINTERFACE_H
#define GRIDGENERATOR_GRIDINTERFACE_H

#include <sfcmm_common.h>
#include "geometry.h"

class GridInterface {
 public:
  GridInterface()          = default;
  virtual ~GridInterface() = default;

  // deleted constructors not needed
  GridInterface(const GridInterface&) = delete;
  GridInterface(GridInterface&&)      = delete;
  auto operator=(const GridInterface&) -> GridInterface& = delete;
  auto operator=(GridInterface&&) -> GridInterface& = delete;

  ////Setter functions.

  /// Set the bounding box for the current grid.
  /// \param bbox Provide the bounding box in the following format {x1, x2, y1, y2...}
  virtual void setBoundingBox(const BoundingBoxInterface& bbox) = 0;

  /// Set the maximum number of cells that this grid can use. Can only be called once!
  /// \param capacity The capacity of this grid object to store cells.
  virtual void setCapacity(GInt capacity) = 0;

  /// Reset the Grid
  virtual void reset() = 0;

  /// Maximum possible level of the grid to store.
  /// \param maxLvl The maximum possible level the grid object can store.
  virtual void setMaxLvl(const GInt maxLvl) = 0;

  virtual void setGeometryManager(std::shared_ptr<GeometryInterface> geom) = 0;

  //// Getter functions.

  /// Get the center of gravity of the grid.
  /// \return Center of gravity of the grid.
  [[nodiscard]] virtual inline auto cog() const -> std::vector<GDouble> = 0;

  /// Length of the sides of the bounding box.
  /// \return Bounding box size.
  [[nodiscard]] virtual inline auto lengthOfBoundingBox() const -> std::vector<GDouble> = 0;

  /// The bounding box of the grid.
  /// \return Bounding box.
  [[nodiscard]] virtual inline auto boundingBox() const -> const BoundingBoxInterface& = 0;

  /// Direction of the largest length of the bounding box.
  /// \return Direction of largest length of the bounding box..
  [[nodiscard]] virtual inline auto largestDir() const -> GInt = 0;

  /// The length of the cell at a given level.
  /// \param lvl The level of the cell length.
  /// \return The length of a cell at the provided level.
  [[nodiscard]] virtual inline auto lengthOnLvl(const GInt lvl) const -> GDouble = 0;

  /// Get the partition level of the grid.
  /// \return Partition level.
  [[nodiscard]] virtual inline auto partitionLvl() const -> GInt = 0;

  /// Set the partition level of the grid.
  /// \return Partition level.
  virtual inline auto partitionLvl() -> GInt& = 0;

  /// The maximum level supported by this grid.
  /// \return The maximum level this grid can possibly have.
  [[nodiscard]] virtual inline auto maxLvl() const -> GInt = 0;

  /// Get the currently highest level present in the grid.
  /// \return Return currently highest level.
  [[nodiscard]] virtual inline auto currentHighestLvl() const -> GInt = 0;

  /// Get the number of dimensions.
  /// \return Dimensions of the grid.
  [[nodiscard]] virtual inline auto dim() const -> GInt = 0;


  //  //// Grid Generation specific
  //  /// Create the grid that is used for partitioning. This grid has the level of the option provided in the grid
  //  /// configuration file. The grid up to this level is always produced on a single MPI rank.
  //  /// \param partioningLvl Level of the partitioing grid.
  //  virtual void createPartitioningGrid(const GInt partioningLvl) = 0;
  //
  //  /// Uniformly refine the grid up to the provided level.
  //  /// \param uniformLvl Level of uniform refinement.
  //  virtual void uniformRefineGrid(const GInt uniformLvl) = 0;
  //
  //  /// Refine the cells that have been marked for refinement.
  //  /// \param noCellsToRefine The number of cells that have been marked.
  //  virtual void refineMarkedCells(const GInt noCellsToRefine) = 0;
  //
  //  /// Mark all boundary cells for refinement
  //  /// \return Number of cells marked for refinement
  //  virtual auto markBndryCells() -> GInt = 0;


  ////IO
  /// Save the grid to a file.
  /// \param fileName The filename(path) of the grid file to be saved.
  /// \param gridOutConfig Json configuration object containing the configuration options of the output format.
  virtual void save(const GString& fileName, const json& gridOutConfig) = 0;

 private:
};
#endif // GRIDGENERATOR_GRIDINTERFACE_H
