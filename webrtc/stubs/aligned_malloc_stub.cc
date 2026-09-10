#include "rtc_base/memory/aligned_malloc.h"
#include <cstdlib>
#include <cstdint>

namespace webrtc {

void* GetRightAlign(const void* ptr, size_t alignment) {
  uintptr_t p = reinterpret_cast<uintptr_t>(ptr);
  return reinterpret_cast<void*>((p + alignment - 1) & ~(alignment - 1));
}

void* AlignedMalloc(size_t size, size_t alignment) {
  void* ptr = nullptr;
  if (posix_memalign(&ptr, alignment, size) != 0)
    return nullptr;
  return ptr;
}

void AlignedFree(void* mem_block) {
  free(mem_block);
}

}  // namespace webrtc
