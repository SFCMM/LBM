// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_SYS_H
#define SFCMM_SYS_H

#include <array>
#include <date.h>
#include <dirent.h>
#include <filesystem>
#include <pwd.h>
#include <unistd.h>

#include "common/sfcmm_types.h"

/// Get the current time as YYYY-MM-DD HH:MM:SS
/// \return Returns the current time stamp as a string.
inline auto dateString() -> GString { return date::format("%F %X", std::chrono::system_clock::now()); }

/// Get the name of the host that the application is running on.
/// \return Returns the name of the host as a string.
inline auto hostString() -> GString {
  static constexpr GInt bufferSize = 1024;

  // Get current hostname
  std::array<char, bufferSize> host_{};
  gethostname(&host_[0], bufferSize - 1);
  host_[bufferSize - 1] = '\0';
  return static_cast<GString>(&host_[0]);
}

/// Get the name of the user that runs the application.
/// \return Returns the name of the user as a string.
inline auto userString() -> GString {
  // Gets the current username
  GString user;

  passwd* p = getpwuid(getuid());
  if(p != nullptr) {
    user = GString(p->pw_name);
  } else {
    user = "n/a";
  }
  return user;
}

/// Get the current working directory.
/// \return Returns the current working directory as a string.
inline auto getCWD() -> GString { return std::filesystem::current_path(); }

/// Check if the given file name is already existing.
/// \param name File name to check for existence.
/// \return File exists?
inline auto isFile(const std::string& name) -> GBool { return std::filesystem::exists(name); }

/// Check if the given path is valid.
/// \param path to check
/// \param generateDir if path doesn't exist create dir (default=false)
/// \return Path exists?
inline auto isPath(const std::string& path, const GBool generateDir = false) -> GBool {
  const GBool existPath = std::filesystem::exists(path);
  if(!existPath && generateDir) {
    return std::filesystem::create_directory(path);
  }
  return existPath;
}

/// Get size of a file in bytes.
/// \param name File name of the file to be sized.
/// \return Number of bytes contained in the file.
inline auto fileSize(const std::string& name) -> GInt { return static_cast<GInt>(std::filesystem::file_size(name)); }
#endif // SFCMM_SYS_H
