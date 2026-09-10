#ifndef RTC_BASE_CHECKS_H_
#define RTC_BASE_CHECKS_H_

#include <cassert>
#include <cstdlib>

namespace rtc {
class NullStream {
 public:
  template <typename T>
  NullStream& operator<<(const T&) { return *this; }
};

struct NullStreamSink {
  void operator&(const NullStream&) {}
};
}  // namespace rtc

#define RTC_DCHECK_IS_ON 0

#define RTC_EAT_STREAM_PARAMS() \
  true ? (void)0 : ::rtc::NullStreamSink() & ::rtc::NullStream()

#define RTC_CHECK(cond) (cond) ? (void)0 : ::rtc::NullStreamSink() & (abort(), ::rtc::NullStream())
#define RTC_CHECK_EQ(a, b) ((a) == (b)) ? (void)0 : ::rtc::NullStreamSink() & (abort(), ::rtc::NullStream())
#define RTC_CHECK_NE(a, b) ((a) != (b)) ? (void)0 : ::rtc::NullStreamSink() & (abort(), ::rtc::NullStream())
#define RTC_CHECK_LE(a, b) ((a) <= (b)) ? (void)0 : ::rtc::NullStreamSink() & (abort(), ::rtc::NullStream())
#define RTC_CHECK_LT(a, b) ((a) < (b)) ? (void)0 : ::rtc::NullStreamSink() & (abort(), ::rtc::NullStream())
#define RTC_CHECK_GE(a, b) ((a) >= (b)) ? (void)0 : ::rtc::NullStreamSink() & (abort(), ::rtc::NullStream())
#define RTC_CHECK_GT(a, b) ((a) > (b)) ? (void)0 : ::rtc::NullStreamSink() & (abort(), ::rtc::NullStream())
#define RTC_CHECK_NOTREACHED() RTC_CHECK(false)

#define RTC_DCHECK(cond) RTC_EAT_STREAM_PARAMS()
#define RTC_DCHECK_EQ(a, b) RTC_EAT_STREAM_PARAMS()
#define RTC_DCHECK_NE(a, b) RTC_EAT_STREAM_PARAMS()
#define RTC_DCHECK_LE(a, b) RTC_EAT_STREAM_PARAMS()
#define RTC_DCHECK_LT(a, b) RTC_EAT_STREAM_PARAMS()
#define RTC_DCHECK_GE(a, b) RTC_EAT_STREAM_PARAMS()
#define RTC_DCHECK_GT(a, b) RTC_EAT_STREAM_PARAMS()
#define RTC_DCHECK_NOTREACHED() RTC_EAT_STREAM_PARAMS()
#define RTC_DCHECK_RUNS_SERIALIZED(x) ((void)0)

#define RTC_NOTREACHED() RTC_CHECK(false)

#endif  // RTC_BASE_CHECKS_H_

#ifndef FRIEND_TEST_ALL_PREFIXES
#define FRIEND_TEST_ALL_PREFIXES(test_case_name, test_name)
#endif
