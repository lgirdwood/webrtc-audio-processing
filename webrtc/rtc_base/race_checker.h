#ifndef RTC_BASE_RACE_CHECKER_H_
#define RTC_BASE_RACE_CHECKER_H_

namespace rtc {
class RaceChecker {
 public:
  constexpr RaceChecker() = default;
};

#define RTC_CHECK_RUNS_SERIALIZED(x) ((void)0)
#define RTC_DCHECK_RUNS_SERIALIZED(x) ((void)0)

}  // namespace rtc

#endif  // RTC_BASE_RACE_CHECKER_H_
