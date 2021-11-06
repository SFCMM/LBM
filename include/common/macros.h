// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_MACROS_H
#define SFCMM_MACROS_H

#ifdef CLANG_COMPILER
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma ide diagnostic   ignored "cppcoreguidelines-macro-usage"
#endif


/// Define macros to stringify literal and expanded macro arguments
///
/// STRINGIFY() can be used to stringify literal macro arguments (e.g.
/// STRINGIFY(__FILE__) becomes "__FILE__"), while XSTRINGIFY() will expand a
/// macro first (e.g. XSTRINGIFY(__FILE__) becomes "macros.h").
#define STRINGIFY(s)  #s
#define XSTRINGIFY(s) STRINGIFY(s)

/// Define a short-hand macros for the location in the code (<file>:<line>)
#define LOC_ __FILE__ ":" XSTRINGIFY(__LINE__)

#define FUN_       static_cast<const char*>(__PRETTY_FUNCTION__)
#define SHORT_FUN_ static_cast<const char*>(__FUNCTION__)

#define AT_ std::string(FUN_) + " (" + LOC_ + ")"

#define __FUNCTION_LOCATION__ std::string(__FILE__) + ": " + std::string(SHORT_FUN_)

#ifndef USE_ASSERTS
#define ASSERT(condition, message)                                                                                                         \
  do {                                                                                                                                     \
  } while(false && (condition))
#endif

#ifdef CLANG_COMPILER
#pragma clang diagnostic pop
#endif
#endif // SFCMM_MACROS_H
