#ifndef LBM_VARIABLES_H
#define LBM_VARIABLES_H

/// Class to keep the positions of variables consistent
/// \tparam EQ The equation that is in use.
/// \tparam NDIM The number of dimensions that are used.
template <LBEquation EQ, GInt NDIM>
class LBMVariables {
 public:
  /// Get the position of the x-axis velocity
  /// \return Position in the variable array for U
  static constexpr auto U() -> GInt { return 0; }

  /// Get the position of the y-axis velocity
  /// \return Position in the variable array for V
  static constexpr auto V() -> GInt {
    if(NDIM > 1) {
      return 1;
    }
    return -1;
  }

  /// Get the position of the z-axis velocity
  /// \return Position in the variable array for W
  static constexpr auto W() -> GInt {
    if(NDIM > 2) {
      return 2;
    }
    return -1;
  }

  /// Get an array of all the velocities for iterations
  /// \return Array of velocity position in variable array
  static constexpr auto velocities() -> std::array<GInt, NDIM> {
    if constexpr(NDIM == 1) {
      return {U()};
    }
    if constexpr(NDIM == 2) {
      return {U(), V()};
    }
    if constexpr(NDIM == 3) {
      return {U(), V(), W()};
    }
    TERMM(-1, "Not implemented!");
  }

  /// Get the position of the velocity from the axis direction.
  /// \param dir Axis direction (0 = x, 1=y, 2=z...)
  /// \return Position of velocity in axis direction in the variable array.
  static constexpr auto velocity(const GInt dir) -> GInt {
    if(dir == 0) {
      return U();
    }
    if(dir == 1 && NDIM > 1) {
      return V();
    }
    if(dir == 2 && NDIM > 2) {
      return W();
    }
    TERMM(-1, "Not implemented!");
  }

  /// String representation of the velocities
  static constexpr std::array<std::string_view, 3> VELSTR = {"U", "V", "W"};

  /// Get the position of the density.
  /// \return Position of rho within the variable array.
  static constexpr auto rho() -> GInt {
    //    if(EQ == LBEquation::Poisson) {
    //      TERMM(-1, "Invalid variable used");
    //    }
    return velocity(NDIM - 1) + 1;
  }


  // todo: fix for thermal
  static constexpr auto temperature() -> GInt { return rho() + 1; }

  // todo: fix for thermal
  static constexpr auto electricPotential() -> GInt {
    if(EQ == LBEquation::Poisson) {
      return 0;
    }
    return rho() + 1;
  }

  static constexpr auto varStr(const GInt VID) -> std::string_view {
    if(VID == electricPotential()) {
      return "P";
    }

    if(VID == U()) {
      return "U";
    }

    if(VID == V()) {
      return "V";
    }

    if(VID == W()) {
      return "W";
    }

    if(VID == rho()) {
      return "rho";
    }
    TERMM(-1, "Invalid VID");
  }
};


#endif // LBM_VARIABLES_H
