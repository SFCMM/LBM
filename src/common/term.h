#ifndef GRIDGENERATOR_TERM_H
#define GRIDGENERATOR_TERM_H
#include "common/log.h"
#include "common/util/backtrace.h"

[[noreturn]] inline void term(const GInt errorCode, const GString& location, const GString& message = "") {
  if(errorCode != 0) {
    std::stringstream s;
    s << "\n";
    s << "Rank " << MPI::globalDomainId() << " threw exit code " << errorCode << "\n";
    s << "Error in " << location << ": " << message << "\n";
    s << "\n"
      << "Program is aborting!!\n";
    std::cerr << s.str() << std::flush;
    logger << s.str() << std::endl;

    // Print backtrace (if enabled)
    BACKTRACE();

    // Close the log file to make sure that no MPI error occurs from the
    // unclosed file, and that a proper XML footer is written
    logger.close();
    MPI_Abort(MPI_COMM_WORLD, static_cast<int>(errorCode));
  } else {
    // Close the log file to make sure that no MPI error occurs from the
    // unclosed file, and that a proper XML footer is written
    logger.close();
    // Call MPI_Finalize to ensure proper MPI shutdown
    MPI_Finalize();

    // Exit the program
    exit(0);
  }
  exit(static_cast<int>(errorCode));
}

#define TERMM(exitval, msg)                                                                                                                \
  do {                                                                                                                                     \
    term(exitval, AT_, msg);                                                                                                               \
  } while(false)

#ifdef USE_ASSERTS
#define ASSERT(condition, message)                                                                                                         \
  do {                                                                                                                                     \
    if(!(condition)) {                                                                                                                     \
      std::cerr << "Assertion `" #condition "` failed in " << __FILE__ << " line " << __LINE__ << ": " << message << std::endl;            \
      TERMM(1, "ASSERTION FAILED");                                                                                                        \
    }                                                                                                                                      \
  } while(false)
#else
#ifndef ASSERT
#define ASSERT(condition, message)                                                                                                         \
  do {                                                                                                                                     \
  } while(false && (condition))
#endif
#endif

#endif // GRIDGENERATOR_TERM_H
