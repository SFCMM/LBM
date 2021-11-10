#ifndef LBM_CONFIGURATION_H
#define LBM_CONFIGURATION_H
#include <json.h>
using json = nlohmann::json;

class configuration{
 public:

  void setConfiguration(const GString& configFileName){
    m_configFileName = configFileName;
  }

  void setConfiguration(const json& config){
    m_config = config;
    setupUnusedTracking();
  }



  void loadConfiguration(const GString& section=""){
    if(!m_config.empty()){
      logger << "Warning loading configuration, but configuration already loaded!" << std::endl;
      TERMM(-1, "Invalid operation!");
    }
    logger << "Loading configuration file [" << configFile() << "]" << std::endl;

    // 1. open configuration file on root process
    if(MPI::isRoot() && isFile(configFile())) {
      std::ifstream configFileStream(configFile());
      configFileStream >> m_config;
      if(!section.empty()){
        m_config = m_config[section];
      }
      setupUnusedTracking();
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

  auto has_any_key_value(const GString& key, const GString& value) const -> GBool {
    auto has_config_value_key = [=](const json& cc, const GString& key, const GString& val){
      // entry exists
      if(cc.contains(key)){
        // value matches
        if(cc[key] == val){
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
      for(const auto& item : tmp){
        if(item.is_object()) {
          stack.emplace(item);
        }
      }
    } while(!stack.empty());

    return false;
  }

  void unusedConfigValues() {
    GInt i = 0;
    logger << "The following values in the configuration file are unused: \n";
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

  auto config() const -> const json&{
    return m_config;
  }

 private:
  void setupUnusedTracking(){
    // put all available keys in map to keep track of usage
    for(const auto& element : m_config.items()) {
      m_configKeys.emplace(element.key(), false);
    }
  }

  GString                            m_configFileName = "grid.json";
  json                               m_config{};
  std::unordered_map<GString, GBool> m_configKeys{};
};

#endif // LBM_CONFIGURATION_H
