// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_TIMER_H
#define SFCMM_TIMER_H

#include <chrono>
#include <cmath>
#include <common/math/mathfunctions.h>
#include <ctime>
#include <iomanip>
#include <string>
#include <sys/times.h>
#include <unistd.h>
#include <utility>
#include <vector>
#include "sfcmm_types.h"
#include "macros.h"
#ifdef USE_ASSERTS
#include "term.h"
#endif

// todo: documentation
// todo: tests

using chronoTimePoint = std::chrono::time_point<std::chrono::system_clock, std::chrono::duration<double>>;

class Timer {
 public:
  Timer(std::string n, const GInt g, const GInt id, const GInt p)
    : m_name(std::move(n)),
      m_group(g),
      m_timerId(id),
      m_parent(p),
      m_recordedTime(0),
      m_status(Timer::Uninitialized),
      m_subTimers(0),
      m_display(false) {}
  enum { Uninitialized = 0, Running = 1, Stopped = 0 };

  void start(chronoTimePoint t) {
    m_oldCpuTime = m_cpuTime;
    m_cpuTime    = t;
    m_status     = Timer::Running;
  }

  void stop() { m_status = Timer::Stopped; }

  [[nodiscard]] auto inline name() const -> GString { return m_name; }

  [[nodiscard]] auto inline status() const -> GInt { return m_status; }

  [[nodiscard]] auto inline group() const -> GInt { return m_group; }

  auto inline recordedTime() -> GDouble& { return m_recordedTime; }
  [[nodiscard]] auto inline recordedTime() const -> GDouble { return m_recordedTime; }

  auto inline subTimers() -> std::vector<GInt>& { return m_subTimers; }
  [[nodiscard]] auto inline subTimers() const -> const std::vector<GInt>& { return m_subTimers; }

  [[nodiscard]] auto inline oldCpuTime() const -> const chronoTimePoint& { return m_oldCpuTime; }

  auto inline cpuTime() -> chronoTimePoint& { return m_cpuTime; }
  [[nodiscard]] auto inline cpuTime() const -> const chronoTimePoint& { return m_cpuTime; }

  [[nodiscard]] auto inline display() const -> GBool { return m_display; }
  auto inline display() -> GBool& { return m_display; }

 private:
  GString           m_name;  ///< Timer Name
  GInt              m_group; ///< Group Id
  GInt              m_timerId = -1;
  GInt              m_parent  = -1; ///< Parent timer id
  chronoTimePoint   m_cpuTime;      ///< CPU time
  chronoTimePoint   m_oldCpuTime;   ///< Old CPU time (for timer restart)
  GDouble           m_recordedTime; ///< Time recorded on the timer.
  GInt              m_status;       ///< Timer's status, see enum:
  std::vector<GInt> m_subTimers{};
  GBool             m_display;
};

/// \brief TimerManager manages all Timers and allows primitive profiling.
///
/// Usage:
/// - NEW_TIMER(string name) creates a new timer with name "name" and
/// returns its index (a GInt that you can use to access it).
/// - RESET_TIMER(GInt timerId): resets timerId.
/// - START_TIMER(timerId)/STOP_TIMER(timerId) work as expected.
/// - DISPLAY_TIMER(timerId) writes the timerId name and time to the log.
/// - DISPLAY_ALL_TIMERS: displays all timers.
///
/// Example:
///
/// RESET_TIMER;
/// GInt myTimer;
/// myTimer = NEW_TIMER("My timer");
///
/// START_TIMER(myTimer);
/// f1(); // Function will be timed.
/// STOP_TIMER(myTimer);
/// f2(); // Function will not be timed.
/// START_TIMER(myTimer);
/// f3(); // Function will be timed.
/// STOP_TIMER(myTimer);
///
/// DISPLAY_TIMER(myTimer);
///
///
class TimerManager {
  friend auto timers() -> TimerManager& {
    static TimerManager timers;
    return timers;
  }

 public:
  inline auto newGroup(const std::string& groupName) -> GInt;
  inline auto newGroupTimer(const std::string& timerName, const GInt groupId) -> GInt;
  inline auto newSubTimer(const std::string& timerName, const GInt timerId) -> GInt;
  inline auto returnTimer(const GInt timerId) -> GDouble;
  inline auto returnTimerTime(const GInt timerId) -> GDouble;
  inline void startTimer(const GInt timerId, const GString& pos = "");       // start
  inline void stopTimer(const GInt timerId, const GString& pos = "");        // stop
  inline void resetTimer(const GInt timerId);                                // reset
  inline void recordTimer(const GInt timerId);                               // record
  inline void recordTimerStart(const GInt timerId, const GString& pos = ""); // reset + start
  inline void recordTimerStop(const GInt timerId, const GString& pos = "");  // stop + record
  inline void recordTimers();
  inline void stopAllTimers();
  inline void stopAllRecordTimers(); // stop all + record stopped
  inline void displayTimer(const GInt timerId);
  inline void displayTimerNoToggleDisplayed(const GInt timerId);
  inline void displayAllTimers();
  inline void displayAllTimers(const GInt groupId);
  inline void resetTimers();
  inline void resetRecord(const GInt timerId);
  inline void resetRecords();

  // delete: copy construction, and copy assignment
  TimerManager(TimerManager&) = delete;
  auto operator=(const TimerManager&) -> TimerManager& = delete;
  auto operator=(TimerManager&&) -> TimerManager& = delete;
  TimerManager(const TimerManager&)               = delete;
  TimerManager(TimerManager&&)                    = delete;

 private:
  TimerManager()  = default;
  ~TimerManager() = default;

  std::vector<std::string> m_groups;
  std::vector<Timer>       m_timers;

  static inline auto time() -> chronoTimePoint;
  inline void displayTimer_(const GInt timerId, const GBool toggleDisplayed = true, const GInt tIndent = 0, const GDouble superTime = -1.0);
  inline void displayTimerHeader_();
  inline void displayTimerGroupHeader_(const GInt groupId);
  [[nodiscard]] static inline auto indent(const GInt pIndent) -> GInt { return pIndent + 2; };
};

auto timers() -> TimerManager&;

inline auto TimerManager::time() -> chronoTimePoint {
  using clock = std::chrono::system_clock;
  return clock::now();
}

inline void TimerManager::resetRecord(const GInt timerId) { m_timers[timerId].recordedTime() = 0.0; }

inline void TimerManager::resetRecords() {
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    resetRecord(static_cast<GInt>(timerId));
  }
}

inline void TimerManager::recordTimerStop(const GInt timerId, const GString& pos) {
  if(timerId < 0) {
    return;
  }
#ifdef TIMER_SYNC
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  stopTimer(timerId, pos);
  recordTimer(timerId);
}

inline void TimerManager::recordTimer(const GInt timerId) { m_timers[timerId].recordedTime() += returnTimer(timerId); }
inline void TimerManager::recordTimers() {
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    recordTimer(static_cast<GInt>(timerId));
  }
}

inline void TimerManager::resetTimer(const GInt timerId) { m_timers[timerId].stop(); }

inline void TimerManager::resetTimers() {
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    resetTimer(static_cast<GInt>(timerId));
  }
}

/// Creates a new timer group and returns its groupId.
inline auto TimerManager::newGroup(const std::string& name) -> GInt {
  m_groups.push_back(name);
  return static_cast<GInt>(m_groups.size() - 1);
}

/// Creates a new timer and returns its timerId.
inline auto TimerManager::newGroupTimer(const std::string& name, const GInt groupId) -> GInt {
  ASSERT(static_cast<std::size_t>(groupId) < m_groups.size() && groupId > -1,
         "groupId: " + std::to_string(groupId) + " does not exists | name: " + name);
  const GInt newTimerId = static_cast<GInt>(m_timers.size());
  m_timers.emplace_back(name, groupId, newTimerId, -1);
  return newTimerId;
}

/// Creates a new timer and returns its timerId.
inline auto TimerManager::newSubTimer(const std::string& name, const GInt timerId) -> GInt {
  if(timerId < 0) {
    return -1;
  }

  ASSERT(static_cast<std::size_t>(timerId) < m_timers.size(),
         "timerId " + std::to_string(timerId) + " does not exist when trying to create subtimer with name " + name);

  const GInt groupId    = m_timers[timerId].group();
  const GInt newTimerId = static_cast<GInt>(m_timers.size());
  m_timers.emplace_back(name, groupId, newTimerId, timerId);
  m_timers[timerId].subTimers().push_back(newTimerId);
  return newTimerId;
}

inline void TimerManager::recordTimerStart(const GInt timerId, const GString& pos) {
  if(timerId < 0) {
    return;
  }
#ifdef TIMER_SYNC
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  resetTimer(timerId);
  startTimer(timerId, pos);
}

inline void TimerManager::startTimer(const GInt timerId, const GString& pos) {
  ASSERT(m_timers[timerId].status() != Timer::Running,
         "The timer " + m_timers[timerId].name() + " with id: " + std::to_string(timerId)
             + " can't be started because it is already running! " + pos);
  std::ignore = pos;

  chronoTimePoint t = time();
  m_timers[timerId].start(t);
}

/// Returns the timer Value.
inline auto TimerManager::returnTimer(const GInt timerId) -> GDouble {
#ifdef TIMER_RANK_AVG
  const GDouble t       = m_timers[timerId].cpuTime().time_since_epoch().count();
  GDouble       tmp_rcv = 0.0;
  MPI_Reduce(&t, &tmp_rcv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  return tmp_rcv / static_cast<GDouble>(MPI::globalNoDomains());
#else
  return m_timers[timerId].cpuTime().time_since_epoch().count();
#endif
}

// Returns the recorded time
inline auto TimerManager::returnTimerTime(const GInt timerId) -> GDouble { return m_timers[timerId].recordedTime(); }

/// Stops the timer and sets its final value.
inline void TimerManager::stopTimer(const GInt timerId, const GString& pos) {
  if(timerId < 0) {
    return;
  }

  if(m_timers[timerId].status() == Timer::Running) {
    const chronoTimePoint t = time();

    m_timers[timerId].cpuTime() = chronoTimePoint(t - m_timers[timerId].cpuTime());

    m_timers[timerId].stop();

    std::ignore = pos;
  } else {
    const GString msg = "The timer '" + m_timers[timerId].name() + "' can't be stopped because it is not running! " + pos;
    std::cerr << msg << std::endl;
  }
}

// Stops all timers.
inline void TimerManager::stopAllTimers() {
  for(std::size_t i = 0, e = m_timers.size(); i != e; ++i) {
    if(m_timers[i].status() == Timer::Running) {
      stopTimer(static_cast<GInt>(i), AT_);
    }
  }
}

// Stops all timers and record the timers that were stopped
inline void TimerManager::stopAllRecordTimers() {
  // for(std::size_t i = 0, e = m_timers.size(); i != e; ++i) {
  for(GInt i = static_cast<GInt>(m_timers.size() - 1); i >= 0; i--) {
    if(m_timers[i].status() == Timer::Running) {
#ifdef TIMER_SYNC
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      stopTimer(i, AT_);
      recordTimer(i);
    }
  }
}

inline void TimerManager::displayTimer_(const GInt timerId, const GBool toggleDisplayed, const GInt tIndent, const GDouble superTime) {
  GBool running = false;
  if(m_timers[timerId].display()) {
    return;
  }

  if(m_timers[timerId].status() == Timer::Running) {
    running = true;
    stopTimer(timerId, AT_);
  }
  logger.width(50);
  logger.setf(std::ios::left);
  std::stringstream indentedName;

  // Calculate time relative to the parent timer
  GDouble percentage = NAN;
  if(superTime < 0.0) {
    // If the parent time is less than zero, that means that there is no parent
    // timer and the percentage should be 100%
    percentage = 100.0;
  } else if(approx(superTime, 0.0, GDoubleEps)) {
    // If the parent time is approximately zero, that probably means that the
    // timer was never run - therefore the percentage is set to 0%
    percentage = 0.0;
  } else {
    // Otherwise calculate the percentage as the fraction of this timer vs. the
    // parent timer times 100%
    percentage = 100.0 * m_timers[timerId].recordedTime() / superTime;
  }

  indentedName << std::string(tIndent, ' ');
  indentedName << "[" << std::fixed << std::setprecision(1) << std::setw(4) << std::setfill('0') << std::right << percentage << std::left
               << "%] ";
  indentedName << m_timers[timerId].name();
  logger << indentedName.str() << std::right;
  logger.precision(6);
  logger.width(20);
  logger << m_timers[timerId].recordedTime() << std::left << " [sec]";
  if(toggleDisplayed) {
    m_timers[timerId].display() = true;
  }
  logger << std::endl;
  for(std::size_t sub = 0, last = m_timers[timerId].subTimers().size(); sub < last; ++sub) {
    const GInt new_indent = indent(tIndent);
    displayTimer_(m_timers[timerId].subTimers()[sub], toggleDisplayed, new_indent, m_timers[timerId].recordedTime());
  }
  if(running) {
    startTimer(timerId);
  }
}

inline void TimerManager::displayTimerHeader_() {}

inline void TimerManager::displayTimerGroupHeader_(const GInt groupId) {
  logger << "------------------------------------------------------------------"
            "--------------"
         << std::endl;
  logger.width(50);
  logger.precision(12);
  logger.setf(std::ios::left);
  logger << "Group";
  logger.width(40);
  logger << m_groups[groupId] << std::endl;
}

inline void TimerManager::displayAllTimers() {
  ASSERT(!m_timers.empty(), "ERROR: no timers have been created!");
  for(std::size_t groupId = 0, e = m_groups.size(); groupId != e; ++groupId) {
    displayAllTimers(static_cast<GInt>(groupId));
  }
}

inline void TimerManager::displayAllTimers(const GInt groupId) {
  ASSERT(static_cast<GInt>(m_timers.size()) > 0, "ERROR: no timers have been created!");
  ASSERT(static_cast<std::size_t>(groupId) < m_groups.size() && groupId > -1, "ERROR: groupId does not exists");
  for(auto& m_timer : m_timers) {
    m_timer.display() = false;
  }
  displayTimerGroupHeader_(groupId);
  displayTimerHeader_();
  for(std::size_t timerId = 0, e = m_timers.size(); timerId != e; ++timerId) {
    if(m_timers[timerId].group() == groupId) {
      displayTimer_(static_cast<GInt>(timerId));
    }
  }
  for(auto& m_timer : m_timers) {
    m_timer.display() = false;
  }
}

inline void TimerManager::displayTimer(const GInt timerId) {
  ASSERT(static_cast<std::size_t>(timerId) < m_timers.size(), "ERROR: timer timerId does not exist");
  displayTimerHeader_();
  displayTimer_(timerId);
}

inline void TimerManager::displayTimerNoToggleDisplayed(const GInt timerId) {
  ASSERT(static_cast<std::size_t>(timerId) < m_timers.size(), "ERROR: timer timerId does not exist");
  displayTimerHeader_();
  displayTimer_(timerId, false);
}

inline void NEW_TIMER_GROUP_NOCREATE(GInt& id, const GString& groupName) { id = timers().newGroup(groupName); }

inline void NEW_TIMER_NOCREATE(GInt& id, const GString& timerName, const GInt groupId) { id = timers().newGroupTimer(timerName, groupId); }

inline void NEW_SUB_TIMER_NOCREATE(GInt& id, const GString& timerName, const GInt timerId) {
  id = timers().newSubTimer(timerName, timerId);
}

inline void START_TIMER(const GInt timerId) { timers().startTimer(timerId); }

inline void RETURN_TIMER(const GInt timerId) { timers().returnTimer(timerId); }

inline void RETURN_TIMER_TIME(const GInt timerId) { timers().returnTimerTime(timerId); }

inline void STOP_TIMER(const GInt timerId) { timers().stopTimer(timerId); }

inline void STOP_ALL_TIMERS() { timers().stopAllTimers(); }

inline void RECORD_TIMER(const GInt timerId) { timers().recordTimer(timerId); }

inline void RECORD_TIMERS() { timers().recordTimers(); }

inline void STOP_ALL_RECORD_TIMERS() { timers().stopAllRecordTimers(); }

inline void DISPLAY_TIMER(const GInt timerId) { timers().displayTimer(timerId); }

inline void DISPLAY_TIMER_INTERM(const GInt timerId) { timers().displayTimerNoToggleDisplayed(timerId); }

inline void DISPLAY_ALL_GROUP_TIMERS(const GInt groupId) { timers().displayAllTimers(groupId); }

inline void DISPLAY_ALL_TIMERS() { timers().displayAllTimers(); }

inline void RESET_TIMER(const GInt timerId) { timers().resetTimer(timerId); }

inline void RESET_TIMERS() { timers().resetTimers(); }

inline void RESET_RECORD(const GInt timerId) { timers().resetRecord(timerId); }

inline void RESET_ALL_RECORDS() { timers().resetRecords(); }

// Currently, not used.
//#define NEW_TIMER_GROUP(id, groupName)                 const GInt id = timers().newGroup(groupName)
//#define NEW_TIMER(id, timerName, groupId)              const GInt id = timers().newGroupTimer(timerName, groupId)
//#define NEW_SUB_TIMER(id, timerName, timerId)          const GInt id = timers().newSubTimer(timerName, timerId)
//#define NEW_TIMER_GROUP_STATIC(id, groupName)          static const GInt id = timers().newGroup(groupName)
//#define NEW_TIMER_STATIC(id, timerName, groupId)       static const GInt id = timers().newGroupTimer(timerName, groupId)
//#define NEW_SUB_TIMER_STATIC(id, timerName, timerId)   static const GInt id = timers().newSubTimer(timerName, timerId)

// Macro for debugging purpose (AT_)!
#define RECORD_TIMER_START(timerId) timers().recordTimerStart(timerId, AT_)
#define RECORD_TIMER_STOP(timerId)  timers().recordTimerStop(timerId, AT_)

class TimerProfiling;
class FunctionTiming;

inline auto cpuTime() -> clock_t { return clock(); }

inline auto wallTime() -> GDouble { return MPI_Wtime(); }

class FunctionTiming {
 public:
  explicit FunctionTiming(std::string name)
    : m_initCpuTime(cpuTime()),
      m_deltaCpuTime(0),
      m_tmpCpuTime(0),
      m_initWallTime(wallTime()),
      m_deltaWallTime(0.0),
      m_tmpWallTime(-1.0),
      m_name(std::move(name)){};
  ~FunctionTiming() { m_name = "<deleted>"; };
  auto operator=(const FunctionTiming& t) -> FunctionTiming& = default;
  auto operator=(FunctionTiming&&) -> FunctionTiming& = default;
  FunctionTiming(FunctionTiming&)                     = default;
  FunctionTiming(const FunctionTiming&)               = default;
  FunctionTiming(FunctionTiming&&)                    = default;

  void in() {
    m_tmpCpuTime  = cpuTime();
    m_tmpWallTime = wallTime();
  }

  void out() {
    if(m_tmpCpuTime > 0) {
      m_deltaCpuTime += (cpuTime() - m_tmpCpuTime);
    }
    if(m_tmpWallTime > 0.0) {
      m_deltaWallTime += (wallTime() - m_tmpWallTime);
    }
    m_tmpCpuTime  = -1;
    m_tmpWallTime = -1.0;
  }

  [[nodiscard]] auto getInitCpuTime() const -> clock_t { return m_initCpuTime; }
  [[nodiscard]] auto getDeltaCpuTime() const -> clock_t { return m_deltaCpuTime; }
  [[nodiscard]] auto getInitWallTime() const -> GDouble { return m_initWallTime; }
  [[nodiscard]] auto getDeltaWallTime() const -> GDouble { return m_deltaWallTime; }
  [[nodiscard]] auto getName() const -> std::string { return m_name; }

 private:
  clock_t     m_initCpuTime;
  clock_t     m_deltaCpuTime;
  clock_t     m_tmpCpuTime;
  GDouble     m_initWallTime;
  GDouble     m_deltaWallTime;
  GDouble     m_tmpWallTime;
  std::string m_name;
};

inline auto operator<(const FunctionTiming& a, const FunctionTiming& b) -> GBool {
  if(a.getDeltaCpuTime() == b.getDeltaCpuTime()) {
    return (a.getInitCpuTime() < b.getInitCpuTime());
  }
  return (a.getDeltaCpuTime() < b.getDeltaCpuTime());
}

namespace profiling {
class Tracer;
}
class TimerProfiling {
  friend profiling::Tracer;

 public:
  explicit TimerProfiling(std::string name) : m_initCpuTime(cpuTime()), m_initWallTime(wallTime()), m_name(std::move(name)) {}
  ~TimerProfiling() {
    using namespace std;

    const clock_t exitCpuTime         = cpuTime();
    const GDouble exitWallTime        = wallTime();
    const GDouble thresholdPercentage = 0.5;
    stringstream  sstream;
    sstream << "    CPU      WALL   FUNCTION                    >> profile: '" << m_name << "' <<";
    const string header = sstream.str();
    for(std::size_t i = 0; i < header.size(); i++) {
      logger << "_";
    }
    logger << endl;
    logger << header << endl;
    for(std::size_t i = 0; i < header.size(); i++) {
      logger << "-";
    }
    logger << endl;
    GInt counter    = 0;
    GInt supCounter = 0;
    if(!s_functionTimings.empty()) {
      sort(s_functionTimings.begin(), s_functionTimings.end());
      reverse(s_functionTimings.begin(), s_functionTimings.end());
      for(auto& s_functionTiming : s_functionTimings) {
        if(s_functionTiming.getInitCpuTime() < m_initCpuTime) {
          continue;
        }
        const GDouble relCpuTime =
            100.0 * getCpuTimeSecs(s_functionTiming.getDeltaCpuTime()) / max(1e-15, getCpuTimeSecs(exitCpuTime - m_initCpuTime));
        const GDouble relWallTime = 100.0 * s_functionTiming.getDeltaWallTime() / max(1e-15, (exitWallTime - m_initWallTime));
        if(relCpuTime < thresholdPercentage) {
          supCounter++;
          continue;
        }
        stringstream ss;
        ss << fixed << setprecision(2) << relCpuTime << "%   " << relWallTime << "%   " << s_functionTiming.getName();
        logger << ss.str() << endl;
        counter++;
      }
      if(supCounter > 0) {
        logger << "  .....     .....   (" << supCounter << " shorter timings with CPU<" << thresholdPercentage << "% were suppressed)"
               << endl;
      }
    }
    if(counter == 0) {
      logger << "No timings recorded for timer '" << m_name << "'." << endl;
    }
    for(std::size_t i = 0; i < header.size(); i++) {
      logger << "-";
    }
    logger << endl;
    logger << "Total cpu time:  " << printTime(getCpuTimeSecs(exitCpuTime - m_initCpuTime)) << endl;
    logger << "Total wall time: " << printTime(exitWallTime - m_initWallTime) << endl;
    for(std::size_t i = 0; i < header.size(); i++) {
      logger << "_";
    }
    logger << endl;
  }

  auto operator=(const TimerProfiling& t) -> TimerProfiling& = delete;
  auto operator=(TimerProfiling&&) -> TimerProfiling& = delete;
  TimerProfiling(TimerProfiling&)                     = delete;
  TimerProfiling(const TimerProfiling&)               = delete;
  TimerProfiling(TimerProfiling&&)                    = delete;

  static auto getTimingId(const std::string& name) -> GInt {
    GInt tId = -1;
    if(!s_functionTimings.empty()) {
      for(std::vector<FunctionTiming>::size_type i = 0; i < s_functionTimings.size(); i++) {
        if(s_functionTimings[i].getName() == name) {
          tId = static_cast<GInt>(i);
        }
      }
    }
    if(tId < 0) {
      tId = static_cast<GInt>(s_functionTimings.size());
      s_functionTimings.emplace_back(name);
    }
    ASSERT(tId > -1, "Non-existing timer");
    return tId;
  }

  static auto getCpuTimeSecs(clock_t cput) -> GDouble { return (static_cast<GDouble>(cput) / static_cast<GDouble>(CLOCKS_PER_SEC)); }

  static auto printTime(GDouble secs) -> GString {
    using namespace timeconst;
    std::stringstream time;
    time.str("");
    GDouble rem = secs;
    if(rem > DDAY) {
      const GDouble div = floor(rem / DDAY);
      time << (static_cast<GInt>(div)) << " days, ";
      rem -= div * DDAY;
    }
    if(rem > DHOUR) {
      const GDouble div = floor(rem / DHOUR);
      time << (static_cast<GInt>(div)) << " hours, ";
      rem -= div * DHOUR;
    }
    if(rem > DMINUTE) {
      const GDouble div = floor(rem / DMINUTE);
      time << (static_cast<GInt>(div)) << " mins, ";
      rem -= div * DMINUTE;
    }
    time << rem << " secs";
    return time.str();
  }

 private:
  const clock_t                             m_initCpuTime;
  const GDouble                             m_initWallTime;
  const std::string                         m_name;
  inline static std::vector<FunctionTiming> s_functionTimings;
};

/// Helper class for automatic tracing
namespace profiling {
class Tracer {
 public:
  Tracer(const GString& funloc) : m_timerId(TimerProfiling::getTimingId(funloc)) { TimerProfiling::s_functionTimings[m_timerId].in(); }
  ~Tracer() { TimerProfiling::s_functionTimings[m_timerId].out(); }
  // delete: copy construction, and copy assignment
  Tracer(Tracer&) = delete;
  auto operator=(const Tracer&) -> Tracer& = delete;
  auto operator=(Tracer&&) -> Tracer& = delete;
  Tracer(const Tracer&)               = delete;
  Tracer(Tracer&&)                    = delete;

 private:
  GInt m_timerId;
};
} // namespace profiling

#define PROFILE() profiling::Tracer USE_ONLY_ONCE_PER_SCOPE(__FUNCTION_LOCATION__)

#endif // SFCMM_TIMER_H
