// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_BACKTRACE_H
#define SFCMM_BACKTRACE_H

#include <cmath>
#include <iostream>
#include <sstream>
#include "common/compiler_config.h"

// todo: activate

#ifdef GCC_COMPILER
// Needed for stack trace
#include <cxxabi.h>
#include <execinfo.h>
#endif

// General backtrace macro
#if defined(ENABLE_BACKTRACE)
#include "llvm/Support/Signals.h"
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif
#include "llvm/Support/raw_ostream.h"

#define BACKTRACE()                                                                                                                        \
  do {                                                                                                                                     \
    debug::backtrace();                                                                                                                    \
  } while(false)
#else
#define BACKTRACE()                                                                                                                        \
  do {                                                                                                                                     \
  } while(false)
#endif


namespace debug {

/**
 * \brief Prints a backtrace of the function call path if possible.
 *
 * \details Uses the LLVMSupport library to print a stacktrace
 *          Note: use BACKTRACE(...) instead of calling this method directly
 *
 */

#if defined(ENABLE_BACKTRACE)
inline void backtrace() {
  llvm::errs() << "Backtrace (line numbers may be too large by 1-3 lines):\n";
  llvm::sys::PrintStackTrace(llvm::errs());
}
#endif

} // namespace debug

#endif // SFCMM_BACKTRACE_H
