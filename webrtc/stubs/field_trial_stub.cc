#include "system_wrappers/include/field_trial.h"
#include <string>

namespace webrtc {
namespace field_trial {
std::string FindFullName(absl::string_view name) {
    return "";
}
void InitFieldTrialsFromString(const char* trials_string) {}
const char* GetFieldTrialString() { return ""; }
bool FieldTrialsStringIsValid(absl::string_view trials_string) { return true; }
std::string MergeFieldTrialsStrings(absl::string_view first, absl::string_view second) { return ""; }
}  // namespace field_trial
}  // namespace webrtc
