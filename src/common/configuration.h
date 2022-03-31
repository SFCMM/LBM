#ifndef LBM_CONFIGURATION_H
#define LBM_CONFIGURATION_H
#include "json.h"
#include <sfcmm_common.h>
#include <utility>
#include <stack>
#include "term.h"
#include "mathexpr.h"

using json = nlohmann::json;

template <GInt NDIM>
using Point = VectorD<NDIM>;

namespace config {
template <typename T>
static constexpr inline auto required_config_value(const json& config, const GString& key) -> T {
  // todo: check for types
  if(config.template contains(key)) {
    return static_cast<T>(config[key]);
  }
  cerr0 << "CONFIG:" << std::endl << config << std::endl;
  TERMM(-1, "The required configuration value is missing: " + key);
}

template <GInt NDIM>
static inline auto required_config_value(const json& config, const GString& key) -> Point<NDIM> {
  if(config.template contains(key)) {
    // todo: check size
    // value is just a number
    if(NDIM == 1 && config[key].is_number()) {
      GDouble tmp = config[key];
      return Point<NDIM>(&tmp);
    }
    // array could be also just be an array of size 1
    if(config[key].is_array()) {
      std::vector<GDouble> tmp = config[key];
      return Point<NDIM>(tmp.data());
    }
    std::stringstream ss;
    ss << config;
    TERMM(-1, "Not a number! config: " + ss.str() + " for the value: " + key);
  }
  cerr0 << config << std::endl;
  TERMM(-1, "The required configuration value is missing: " + key);
}

template <GInt NDIM>
static inline auto required_math_expression(const json& config, const GString& key) -> std::unique_ptr<MathExpression<NDIM>> {
  if(config.template contains(key)) {
    const GString expr = config[key];
    return std::make_unique<MathExpression<NDIM>>(expr);
  }
  cerr0 << config << std::endl;
  TERMM(-1, "The required configuration value is missing: " + key);
}

template <typename T>
static constexpr inline auto opt_config_value(const json& config, const GString& key, const T& defaultValue) -> T {
  // todo: check for types
  if(config.template contains(key)) {
    return static_cast<T>(config[key]);
  }
  return defaultValue;
}

static inline auto has_config_value(const json& config, const GString& key) -> GBool { return config.contains(key); }
} // namespace config

class ConfigurationAccess;

class Configuration {
 public:
  /// Set the configuration file name
  /// \param configFileName Configuration file name
  void setConfiguration(const GString& configFileName) { m_configFileName = configFileName; }

  /// Set the configuration from json-object
  /// \param config json-object to be used for configuration
  void setConfiguration(const json& config) {
    m_config = config;
    setupUnusedTracking();
  }

  /// Load the configuration from a configuration file based on the object
  /// \param section Section of the configuration file to use (optional default=everything)
  void load(const GString& section = "") {
    if(!m_config.empty()) {
      logger << "Warning loading configuration, but configuration already loaded!" << std::endl;
      TERMM(-1, "Invalid operation!");
    }
    logger << "Loading configuration file [" << configFile() << "]" << std::endl;

    // 1. open configuration file on root process
    if(MPI::isRoot() && isFile(configFile())) {
      std::ifstream configFileStream(configFile());
      configFileStream >> m_config;
      if(!section.empty()) {
        m_config = m_config[section];
      }
      setupUnusedTracking();
      configFileStream.close();
    }
  }

  /// Get a required configuration value (exit if it doesn't exist)
  /// \tparam T Type of the value
  /// \param key Key of the value
  /// \return Configuration value if it exist or exit
  template <typename T>
  [[nodiscard]] auto required_config_value(const GString& key) -> T {
    // todo: check for types
    if(m_config.template contains(key)) {
      m_unusedKeys[key] = true;
      return static_cast<T>(m_config[key]);
    }
    TERMM(-1, "The required configuration value is missing: " + key);
  }

  /// Get a required configuration value (exit if it doesn't exist) [const version]
  /// \tparam T Type of the value
  /// \param key Key of the value
  /// \return Configuration value if it exist or exit
  template <typename T>
  [[nodiscard]] auto required_config_value(const GString& key) const -> T {
    // todo: check for types
    if(m_config.template contains(key)) {
      return static_cast<T>(m_config[key]);
    }
    TERMM(-1, "The required configuration value is missing: " + key);
  }

  /// Get a required configuration value (exit if it doesn't exist) [Point Version]
  /// \tparam NDIM Dimensionality of the point
  /// \param key Key of the value
  /// \return Configuration value if it exist or exit
  template <GInt NDIM>
  [[nodiscard]] auto required_config_value(const GString& key) -> Point<NDIM> {
    // todo: check for types
    if(m_config.template contains(key)) {
      m_unusedKeys[key] = true;
      // todo: check size
      std::vector<GDouble> tmp = m_config[key];
      return Point<NDIM>(tmp.data());
    }
    cerr0 << m_config << std::endl;
    TERMM(-1, "The required configuration value is missing: " + key);
  }

  /// Get a required configuration value (exit if it doesn't exist)
  /// \tparam T Type of the value
  /// \param key Key of the value
  /// \param parentObj Parent object path
  /// \return Configuration value if it exist or exit
  template <typename T>
  [[nodiscard]] auto required_config_value(const std::vector<GString>& parentObjPath, const GString& key) -> T {
    // todo: check for types
    json conf = m_config;
    for(const auto& parentKey : parentObjPath) {
      conf = conf[parentKey];
    }
    if(conf.template contains(key)) {
      m_unusedKeys[key] = true;
      return static_cast<T>(conf[key]);
    }
    cerr0 << conf << std::endl;
    TERMM(-1, "The required configuration value is missing: " + key);
  }

  /// Get an optional configuration value (return default value if value not found)
  /// \tparam T Type of the value
  /// \param key Key of the value
  /// \param defaultValue Default value
  /// \return Configuration value if it exists or the given default
  template <typename T>
  [[nodiscard]] auto opt_config_value(const GString& key, const T& defaultValue) -> T {
    // todo: check for types
    if(m_config.template contains(key)) {
      m_unusedKeys[key] = true;
      return static_cast<T>(m_config[key]);
    }
    return defaultValue;
  }

  /// Check if a configuration value exists
  /// \param key Key of the value to check
  /// \return true if it exists
  auto has_config_value(const GString& key) -> GBool { return m_config.template contains(key); }

  /// Exists the provided key-value pair in the configuration?
  /// \param key The key of the value
  /// \param value The value
  /// \return The key-value pair exits -> true
  auto has_any_key_value(const GString& key, const GString& value) const -> GBool {
    auto has_config_value_key = [=](const json& configurationObj, const GString& _key, const GString& val) {
      // entry exists
      if(configurationObj.contains(_key)) {
        // value matches
        if(configurationObj[_key] == val) {
          return true;
        }
      }
      return false;
    };

    std::stack<json> stack;
    stack.emplace(m_config);

    do {
      json tmp = stack.top();
      stack.pop();

      if(has_config_value_key(tmp, key, value)) {
        return true;
      }
      for(const auto& item : tmp) {
        if(item.is_object()) {
          stack.emplace(item);
        }
      }
    } while(!stack.empty());

    return false;
  }

  /// Return the json configuration object of a certain name and mark it as used
  /// \param object Object to get
  /// \param pObject Parent object
  /// \return Object with the name [object]
  auto getObject(const GString& object, const GString& pObject = "") -> json {
    const json& conf = (pObject.empty()) ? m_config : m_config[pObject];
    if(!pObject.empty()) {
      m_unusedKeys[pObject] = true;
    }
    return conf[object];
  }

  /// Get all objects that are within the parentobject
  /// \param pObject Parent object name.
  /// \return List of all the object identifiers.
  auto getAllObjects(const GString& pObject = "") const -> std::vector<GString> {
    std::vector<GString> tmp_objL;
    const json&          conf = (pObject.empty()) ? m_config : m_config[pObject];
    for(const auto& [key, item] : conf.items()) {
      if(item.is_object()) {
        tmp_objL.emplace_back(key);
      }
    }
    return tmp_objL;
  }

  /// Get all the keys with a certain value
  /// \param value Value to search for
  /// \return List of keys with the value
  auto get_all_items_with_value(const GString& value) const -> std::vector<json> {
    auto has_config_value = [=](const json& configObj, const GString& val) {
      // entry exists
      return std::any_of(configObj.begin(), configObj.end(), [&](const auto& item) { return item == val; });
    };

    std::vector<json>   found;
    std::stack<json>    stack;
    std::stack<GString> keyStack;
    stack.emplace(m_config);
    keyStack.emplace("default");


    do {
      json tmp = stack.top();
      stack.pop();
      GString tmpK = keyStack.top();
      keyStack.pop();

      if(has_config_value(tmp, value)) {
        json hit;
        hit[tmpK] = tmp;
        found.emplace_back(hit);
        continue;
      }
      for(const auto& [key, v] : tmp.items()) {
        if(v.is_object()) {
          stack.emplace(v);
          keyStack.emplace(key);
        }
      }
    } while(!stack.empty());

    return found;
  }

  /// Print out a list of unused configuration values
  void unusedConfigValues() {
    GBool none                 = true;
    GInt  noUnusedConfigValues = 0;
    logger << "The following values in the configuration file are unused: \n";
    for(const auto& configKey : m_unusedKeys) {
      if(!configKey.second) {
        logger << "[" << ++noUnusedConfigValues << "] " << configKey.first << "\n";
        none = false;
      }
    }
    if(none) {
      logger << "None \n";
    }
    logger << std::endl;
  }

  /// Get the configuration file name
  /// \return Configuration file name
  auto configFile() const -> GString { return m_configFileName; }

  /// Get the json-object of the configuration
  /// \return JSON-object
  auto config() const -> const json& { return m_config; }

  /// Get an accessor to the configuration limited to the subobject.
  /// \param subobject The subobject to be given access to
  /// \return ConfigurationAccess object
  auto getAccessor(const GString& subobject) -> std::shared_ptr<ConfigurationAccess> {
    if(has_config_value(subobject)) {
      m_bondedConfAccess.emplace_back(std::make_shared<ConfigurationAccess>(subobject, this));
      return m_bondedConfAccess.back();
    }
    return nullptr;
  }

 private:
  /// Init the tracking of unused configuration values
  void setupUnusedTracking() {
    // put all available keys in map to keep track of usage
    for(const auto& element : m_config.items()) {
      m_unusedKeys.emplace(element.key(), false);
    }
  }

  // file name of the configuration
  GString m_configFileName = "grid.json";
  // json configuration object
  json m_config{};
  // map used for tracking the unused values
  std::unordered_map<GString, GBool> m_unusedKeys{};
  // bonded access objects
  std::vector<std::shared_ptr<ConfigurationAccess>> m_bondedConfAccess{};
};


class ConfigurationAccess {
 public:
  ConfigurationAccess(GString prefix, Configuration* parentConf) : m_prefix(std::move(prefix)), m_parentConf(parentConf) {}

  [[nodiscard]] auto getAllObjects() const -> std::vector<GString> { return m_parentConf->getAllObjects(m_prefix); }

  /// Get a required configuration value (exit if it doesn't exist)
  /// \tparam T Type of the value
  /// \param key Key of the value
  /// \param parentObjKey Parent object key string
  /// \return Configuration value if it exist or exit
  template <typename T>
  auto required_config_value(const GString& parentObjKey, const GString& key) -> T {
    // todo: check for types
    if(parentObjKey.empty()) {
      TERMM(-1, "Invalid parentObjKey");
    }
    std::vector<GString> access{m_prefix, parentObjKey};
    return m_parentConf->template required_config_value<T>(access, key);
  }

  /// Get a required configuration value (exit if it doesn't exist)
  /// \tparam T Type of the value
  /// \param key Key of the value
  /// \return Configuration value if it exist or exit
  template <typename T>
  auto required_config_value(const GString& key) -> T {
    // todo: check for types
    std::vector<GString> access{m_prefix};
    return m_parentConf->template required_config_value<T>(access, key);
  }

  auto getObject(const GString& object) -> json { return m_parentConf->getObject(object, m_prefix); }

 private:
  const GString  m_prefix;
  Configuration* m_parentConf = nullptr;
};

#endif // LBM_CONFIGURATION_H
