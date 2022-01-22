#ifndef GRIDGENERATOR_GRIDCELL_PROPERTIES_H
#define GRIDGENERATOR_GRIDCELL_PROPERTIES_H
// todo: merge both of these enum classes since the bits don't matter!

/// Cell properties. (Until 63 bits are used this does not increase memory!)
// bitsets are initialized to 0!
enum class CellProperties {
  periodic,              /// cell is a periodic cell
  halo,                  /// cell is a halo cell (a cell only present for communication purpose is a window cell on a neighbouring domain)
  window,                /// cell is a window cell (a cell used for communication with an associated halo cell on a neighbouring domain)
  markedForDeletion,     /// cell is to be deleted
  bndry,                 /// cell is on a boundary
  inside,                /// cell is on the inside of the geometry
  toRefine,              /// cell is to be refined
  partitionLevelShifted, /// cell is on partition level but is not a partition cell!
  partitionCell,         /// cell is a partition cell
  newCell,               /// cell was recently created
  hasBeenCoarsened,      /// cell has been coarsened
  hasBeenRefined,        /// cell has been refined
  marked,                /// cell has been marked
  ghost,
  leaf,
  solid,
  // <<< add new properties here
  NumProperties
};

namespace grid::cell {
/// Converts property name to underlying integer value
constexpr auto p(const CellProperties property) -> std::underlying_type<CellProperties>::type {
  return static_cast<std::underlying_type<CellProperties>::type>(property);
}

using BitsetType = std::bitset<p(CellProperties::NumProperties)>;
} // namespace grid::cell

#endif // GRIDGENERATOR_GRIDCELL_PROPERTIES_H
