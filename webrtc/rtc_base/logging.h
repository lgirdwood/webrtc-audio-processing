#ifndef RTC_BASE_LOGGING_H_
#define RTC_BASE_LOGGING_H_

#include "rtc_base/checks.h"

namespace rtc {
enum LoggingSeverity { LS_VERBOSE, LS_INFO, LS_WARNING, LS_ERROR, LS_NONE };
}  // namespace rtc

using rtc::LoggingSeverity;
using rtc::LS_VERBOSE;
using rtc::LS_INFO;
using rtc::LS_WARNING;
using rtc::LS_ERROR;
using rtc::LS_NONE;

#define RTC_LOG(sev) ::rtc::NullStream()
#define RTC_LOG_V(sev) ::rtc::NullStream()

#endif  // RTC_BASE_LOGGING_H_
