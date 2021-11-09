#ifndef LBM_CONFIGURATION_H
#define LBM_CONFIGURATION_H
#include <json.h>
using json = nlohmann::json;

class configuration{
 public:

  void setConfiguration(const GString& configFileName){
    m_configFileName = configFileName;
  }

  void loadConfiguration(const GString& section=""){
    logger << "Loading configuration file [" << configFile() << "]" << std::endl;

    // 1. open configuration file on root process
    if(MPI::isRoot() && isFile(configFile())) {
      std::ifstream configFileStream(configFile());
      configFileStream >> m_config;
      if(!section.empty()){
        m_config = m_config[section];
      }
      // put all available keys in map to keep track of usage
      for(const auto& element : m_config.items()) {
        m_configKeys.emplace(element.key(), false);
      }
      configFileStream.close();
    }
  }

  template <typename T>
  auto required_config_value(const GString& key) -> T {
    // todo: check for types
    if(m_config.template contains(key)) {
      m_configKeys[key] = true;
      return static_cast<T>(m_config[key]);
    }
    TERMM(-1, "The required configuration value is missing: " + key);
  }

  template <typename T>
  auto opt_config_value(const GString& key, const T& defaultValue) -> T {
    // todo: check for types
    if(m_config.template contains(key)) {
      m_configKeys[key] = true;
      return static_cast<T>(m_config[key]);
    }
    return defaultValue;
  }

  auto has_config_value(const GString& key) -> GBool {
    return m_config.template contains(key);
  }


  void unusedConfigValues() {
    GInt i = 0;
    logger << "The following values in the configuration file are unused:" << std::endl;
    for(const auto& configKey : m_configKeys) {
      if(!configKey.second) {
        logger << "[" << ++i << "] " << configKey.first << "\n";
      }
    }
    logger << std::endl;
  }

  auto configFile() const -> GString {
    return m_configFileName;
  }

 private:
  GString                            m_configFileName = "grid.json";
  json                               m_config{};
  std::unordered_map<GString, GBool> m_configKeys{};
};

#endif // LBM_CONFIGURATION_H
