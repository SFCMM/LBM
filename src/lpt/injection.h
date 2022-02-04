#ifndef LBM_INJECTION_H
#define LBM_INJECTION_H

#include <queue>

template <GInt NDIM>
class Injector : private Configuration {
 public:
  Injector(const json& config) : m_prng(123342896) {
    setConfiguration(config);

    m_location    = required_config_value<NDIM>("injector_location");
    m_orientation = required_config_value<NDIM>("injector_orientation");

    m_openingAngle      = required_config_value<GDouble>("injector_opening_angle");
    m_holeDiameter      = required_config_value<GDouble>("injector_hole_diameter");
    m_massRate          = required_config_value<GDouble>("mass_rate");
    m_injectedDensity   = required_config_value<GDouble>("injected_density");
    m_injectionVelocity = required_config_value<GDouble>("injection_velocity");
  };
  ~Injector() = default;

  /// Integrate the mass and
  /// \param dt time step size
  /// \return The number of new particles to be added
  auto timeStep(const GDouble dt) -> GInt {
    m_mass += m_massRate * dt;

    const GDouble dropletV    = 1.0 / 6.0 * M_PI * gcem::pow(m_holeDiameter, 3);
    const GDouble dropletMass = m_injectedDensity * dropletV;

    // no particle size distribution set just divide the mass up
    if(!m_distributedSize) {
      // rounded down number of droplets since we only want complete droplets
      const GInt newDroplets = gcem::floor(m_mass / dropletMass);

      // update mass with the actual injected number of droplets
      m_mass -= newDroplets * dropletMass;
      return newDroplets;
    }

    //<- particle size distribution is used

    GInt newDroplets = 0;

    // only valid for RR
    const GDouble minDropletMass = 0.3 * m_injectedDensity * 1.0 / 6.0 * M_PI * gcem::pow(0.3 * m_holeDiameter, 3);
    GDouble       lastD          = m_injectionDiameter.empty() ? NAN : m_injectionDiameter.front();
    clear();

    while(m_mass > minDropletMass) {
      const GDouble distedD = std::isnan(lastD) ? randomDist_rosinRammler(m_holeDiameter, m_prng) : lastD;
      const GDouble distedV = 1.0 / 6.0 * M_PI * gcem::pow(distedD, 3);
      const GDouble distedM = m_injectedDensity * distedV;

      lastD = NAN;
      if(distedM > m_mass) {
        m_injectionDiameter.emplace(distedD);
        return newDroplets;
      }

      m_injectionDiameter.emplace(distedD);
      m_mass -= distedM;
      ++newDroplets;
    }
    return newDroplets;
  }

  /// Getter for the next droplet diameter that is to be added
  /// \return Droplet diameter of the next droplet to be injected
  [[nodiscard]] auto dropletDiameter() -> GDouble {
    if(m_distributedSize) {
      const GDouble next = m_injectionDiameter.front();
      m_injectionDiameter.pop();
      return next;
    }
    return m_holeDiameter;
  }

  /// Getter for location
  /// \return Location of the injector
  [[nodiscard]] auto location() const -> const VectorD<NDIM>& { return m_location; }

  /// Getter for location
  /// \return Location of the injector
  [[nodiscard]] auto orientation() const -> const VectorD<NDIM>& { return m_orientation; }

  /// Getter for orientation
  /// \return orientation of the injector
  [[nodiscard]] auto openingAngle() const -> GDouble { return m_openingAngle; }

  /// Getter for the hole diameter of the injector
  /// \return the hole diameter of the injector
  [[nodiscard]] auto holeDiameter() const -> GDouble { return m_holeDiameter; }

  /// Getter for the density of the injected droplets
  /// \return density of the injected droplets
  [[nodiscard]] auto injectedDensity() const -> GDouble { return m_injectedDensity; }

  /// Getter for the injection velocity
  /// \return Initial velocity of the droplets
  [[nodiscard]] auto injectionVelocity() const -> GDouble { return m_injectionVelocity; }


 private:
  /// Clear for next time step
  void clear() {
    std::queue<GDouble> empty;
    std::swap(m_injectionDiameter, empty);
  }

  /// Configuration values
  VectorD<NDIM> m_location;
  VectorD<NDIM> m_orientation;

  GDouble m_openingAngle;
  GDouble m_holeDiameter;
  GDouble m_massRate;
  GDouble m_injectedDensity;
  GDouble m_injectionVelocity;

  GDouble m_mass = 0;

  // todo: make settable
  GBool               m_distributedSize = true;
  std::queue<GDouble> m_injectionDiameter;
  randxor             m_prng;
};

#endif // LBM_INJECTION_H
