/*
 * http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html
 * Copyright Takuya OOURA, 1996-2001
 *
 * You may use, copy, modify and distribute this code for any purpose (include
 * commercial use) and without fee. Please refer to this package when you modify
 * this code.
 *
 * Changes by the WebRTC authors:
 *    - Trivial type modifications.
 *    - Minimal code subset to do rdft of length 128.
 *    - Optimizations because of known length.
 *    - Removed the global variables by moving the code in to a class in order
 *      to make it thread safe.
 *
 *  All changes are covered by the WebRTC license and IP grant:
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "common_audio/third_party/ooura/fft_size_128/ooura_fft.h"

#include <stdint.h>
#include <string.h>

#include "common_audio/third_party/ooura/fft_size_128/ooura_fft_tables_common.h"
#include "rtc_base/system/arch.h"
#include "system_wrappers/include/cpu_features_wrapper.h"


namespace webrtc {

namespace {

#if !(defined(MIPS_FPU_LE) || defined(WEBRTC_HAS_NEON))
static void cft1st_128_C(float* a) {
  const int n = 128;
  int j, k1, k2;
  float wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
  float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

  // The processing of the first set of elements was simplified in C to avoid
  // some operations (multiplication by zero or one, addition of two elements
  // multiplied by the same weight, ...).
  x0r = a[0] + a[2];
  x0i = a[1] + a[3];
  x1r = a[0] - a[2];
  x1i = a[1] - a[3];
  x2r = a[4] + a[6];
  x2i = a[5] + a[7];
  x3r = a[4] - a[6];
  x3i = a[5] - a[7];
  a[0] = x0r + x2r;
  a[1] = x0i + x2i;
  a[4] = x0r - x2r;
  a[5] = x0i - x2i;
  a[2] = x1r - x3i;
  a[3] = x1i + x3r;
  a[6] = x1r + x3i;
  a[7] = x1i - x3r;
  wk1r = rdft_w[2];
  x0r = a[8] + a[10];
  x0i = a[9] + a[11];
  x1r = a[8] - a[10];
  x1i = a[9] - a[11];
  x2r = a[12] + a[14];
  x2i = a[13] + a[15];
  x3r = a[12] - a[14];
  x3i = a[13] - a[15];
  a[8] = x0r + x2r;
  a[9] = x0i + x2i;
  a[12] = x2i - x0i;
  a[13] = x0r - x2r;
  x0r = x1r - x3i;
  x0i = x1i + x3r;
  a[10] = wk1r * (x0r - x0i);
  a[11] = wk1r * (x0r + x0i);
  x0r = x3i + x1r;
  x0i = x3r - x1i;
  a[14] = wk1r * (x0i - x0r);
  a[15] = wk1r * (x0i + x0r);
  k1 = 0;
  for (j = 16; j < n; j += 16) {
    k1 += 2;
    k2 = 2 * k1;
    wk2r = rdft_w[k1 + 0];
    wk2i = rdft_w[k1 + 1];
    wk1r = rdft_w[k2 + 0];
    wk1i = rdft_w[k2 + 1];
    wk3r = rdft_wk3ri_first[k1 + 0];
    wk3i = rdft_wk3ri_first[k1 + 1];
    x0r = a[j + 0] + a[j + 2];
    x0i = a[j + 1] + a[j + 3];
    x1r = a[j + 0] - a[j + 2];
    x1i = a[j + 1] - a[j + 3];
    x2r = a[j + 4] + a[j + 6];
    x2i = a[j + 5] + a[j + 7];
    x3r = a[j + 4] - a[j + 6];
    x3i = a[j + 5] - a[j + 7];
    a[j + 0] = x0r + x2r;
    a[j + 1] = x0i + x2i;
    x0r -= x2r;
    x0i -= x2i;
    a[j + 4] = wk2r * x0r - wk2i * x0i;
    a[j + 5] = wk2r * x0i + wk2i * x0r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j + 2] = wk1r * x0r - wk1i * x0i;
    a[j + 3] = wk1r * x0i + wk1i * x0r;
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j + 6] = wk3r * x0r - wk3i * x0i;
    a[j + 7] = wk3r * x0i + wk3i * x0r;
    wk1r = rdft_w[k2 + 2];
    wk1i = rdft_w[k2 + 3];
    wk3r = rdft_wk3ri_second[k1 + 0];
    wk3i = rdft_wk3ri_second[k1 + 1];
    x0r = a[j + 8] + a[j + 10];
    x0i = a[j + 9] + a[j + 11];
    x1r = a[j + 8] - a[j + 10];
    x1i = a[j + 9] - a[j + 11];
    x2r = a[j + 12] + a[j + 14];
    x2i = a[j + 13] + a[j + 15];
    x3r = a[j + 12] - a[j + 14];
    x3i = a[j + 13] - a[j + 15];
    a[j + 8] = x0r + x2r;
    a[j + 9] = x0i + x2i;
    x0r -= x2r;
    x0i -= x2i;
    a[j + 12] = -wk2i * x0r - wk2r * x0i;
    a[j + 13] = -wk2i * x0i + wk2r * x0r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j + 10] = wk1r * x0r - wk1i * x0i;
    a[j + 11] = wk1r * x0i + wk1i * x0r;
    x0r = x1r + x3i;
    x0i = x1i - x3r;
    a[j + 14] = wk3r * x0r - wk3i * x0i;
    a[j + 15] = wk3r * x0i + wk3i * x0r;
  }
}

static void cftmdl_128_C(float* a) {
  const int l = 8;
  const int n = 128;
  const int m = 32;
  int j0, j1, j2, j3, k, k1, k2, m2;
  float wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
  float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

  for (j0 = 0; j0 < l; j0 += 2) {
    j1 = j0 + 8;
    j2 = j0 + 16;
    j3 = j0 + 24;
    x0r = a[j0 + 0] + a[j1 + 0];
    x0i = a[j0 + 1] + a[j1 + 1];
    x1r = a[j0 + 0] - a[j1 + 0];
    x1i = a[j0 + 1] - a[j1 + 1];
    x2r = a[j2 + 0] + a[j3 + 0];
    x2i = a[j2 + 1] + a[j3 + 1];
    x3r = a[j2 + 0] - a[j3 + 0];
    x3i = a[j2 + 1] - a[j3 + 1];
    a[j0 + 0] = x0r + x2r;
    a[j0 + 1] = x0i + x2i;
    a[j2 + 0] = x0r - x2r;
    a[j2 + 1] = x0i - x2i;
    a[j1 + 0] = x1r - x3i;
    a[j1 + 1] = x1i + x3r;
    a[j3 + 0] = x1r + x3i;
    a[j3 + 1] = x1i - x3r;
  }
  wk1r = rdft_w[2];
  for (j0 = m; j0 < l + m; j0 += 2) {
    j1 = j0 + 8;
    j2 = j0 + 16;
    j3 = j0 + 24;
    x0r = a[j0 + 0] + a[j1 + 0];
    x0i = a[j0 + 1] + a[j1 + 1];
    x1r = a[j0 + 0] - a[j1 + 0];
    x1i = a[j0 + 1] - a[j1 + 1];
    x2r = a[j2 + 0] + a[j3 + 0];
    x2i = a[j2 + 1] + a[j3 + 1];
    x3r = a[j2 + 0] - a[j3 + 0];
    x3i = a[j2 + 1] - a[j3 + 1];
    a[j0 + 0] = x0r + x2r;
    a[j0 + 1] = x0i + x2i;
    a[j2 + 0] = x2i - x0i;
    a[j2 + 1] = x0r - x2r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[j1 + 0] = wk1r * (x0r - x0i);
    a[j1 + 1] = wk1r * (x0r + x0i);
    x0r = x3i + x1r;
    x0i = x3r - x1i;
    a[j3 + 0] = wk1r * (x0i - x0r);
    a[j3 + 1] = wk1r * (x0i + x0r);
  }
  k1 = 0;
  m2 = 2 * m;
  for (k = m2; k < n; k += m2) {
    k1 += 2;
    k2 = 2 * k1;
    wk2r = rdft_w[k1 + 0];
    wk2i = rdft_w[k1 + 1];
    wk1r = rdft_w[k2 + 0];
    wk1i = rdft_w[k2 + 1];
    wk3r = rdft_wk3ri_first[k1 + 0];
    wk3i = rdft_wk3ri_first[k1 + 1];
    for (j0 = k; j0 < l + k; j0 += 2) {
      j1 = j0 + 8;
      j2 = j0 + 16;
      j3 = j0 + 24;
      x0r = a[j0 + 0] + a[j1 + 0];
      x0i = a[j0 + 1] + a[j1 + 1];
      x1r = a[j0 + 0] - a[j1 + 0];
      x1i = a[j0 + 1] - a[j1 + 1];
      x2r = a[j2 + 0] + a[j3 + 0];
      x2i = a[j2 + 1] + a[j3 + 1];
      x3r = a[j2 + 0] - a[j3 + 0];
      x3i = a[j2 + 1] - a[j3 + 1];
      a[j0 + 0] = x0r + x2r;
      a[j0 + 1] = x0i + x2i;
      x0r -= x2r;
      x0i -= x2i;
      a[j2 + 0] = wk2r * x0r - wk2i * x0i;
      a[j2 + 1] = wk2r * x0i + wk2i * x0r;
      x0r = x1r - x3i;
      x0i = x1i + x3r;
      a[j1 + 0] = wk1r * x0r - wk1i * x0i;
      a[j1 + 1] = wk1r * x0i + wk1i * x0r;
      x0r = x1r + x3i;
      x0i = x1i - x3r;
      a[j3 + 0] = wk3r * x0r - wk3i * x0i;
      a[j3 + 1] = wk3r * x0i + wk3i * x0r;
    }
    wk1r = rdft_w[k2 + 2];
    wk1i = rdft_w[k2 + 3];
    wk3r = rdft_wk3ri_second[k1 + 0];
    wk3i = rdft_wk3ri_second[k1 + 1];
    for (j0 = k + m; j0 < l + (k + m); j0 += 2) {
      j1 = j0 + 8;
      j2 = j0 + 16;
      j3 = j0 + 24;
      x0r = a[j0 + 0] + a[j1 + 0];
      x0i = a[j0 + 1] + a[j1 + 1];
      x1r = a[j0 + 0] - a[j1 + 0];
      x1i = a[j0 + 1] - a[j1 + 1];
      x2r = a[j2 + 0] + a[j3 + 0];
      x2i = a[j2 + 1] + a[j3 + 1];
      x3r = a[j2 + 0] - a[j3 + 0];
      x3i = a[j2 + 1] - a[j3 + 1];
      a[j0 + 0] = x0r + x2r;
      a[j0 + 1] = x0i + x2i;
      x0r -= x2r;
      x0i -= x2i;
      a[j2 + 0] = -wk2i * x0r - wk2r * x0i;
      a[j2 + 1] = -wk2i * x0i + wk2r * x0r;
      x0r = x1r - x3i;
      x0i = x1i + x3r;
      a[j1 + 0] = wk1r * x0r - wk1i * x0i;
      a[j1 + 1] = wk1r * x0i + wk1i * x0r;
      x0r = x1r + x3i;
      x0i = x1i - x3r;
      a[j3 + 0] = wk3r * x0r - wk3i * x0i;
      a[j3 + 1] = wk3r * x0i + wk3i * x0r;
    }
  }
}

static void rftfsub_128_C(float* a) {
  const float* c = rdft_w + 32;
  int j1, j2, k1, k2;
  float wkr, wki, xr, xi, yr, yi;

  for (j1 = 1, j2 = 2; j2 < 64; j1 += 1, j2 += 2) {
    k2 = 128 - j2;
    k1 = 32 - j1;
    wkr = 0.5f - c[k1];
    wki = c[j1];
    xr = a[j2 + 0] - a[k2 + 0];
    xi = a[j2 + 1] + a[k2 + 1];
    yr = wkr * xr - wki * xi;
    yi = wkr * xi + wki * xr;
    a[j2 + 0] -= yr;
    a[j2 + 1] -= yi;
    a[k2 + 0] += yr;
    a[k2 + 1] -= yi;
  }
}

static void rftbsub_128_C(float* a) {
  const float* c = rdft_w + 32;
  int j1, j2, k1, k2;
  float wkr, wki, xr, xi, yr, yi;

  a[1] = -a[1];
  for (j1 = 1, j2 = 2; j2 < 64; j1 += 1, j2 += 2) {
    k2 = 128 - j2;
    k1 = 32 - j1;
    wkr = 0.5f - c[k1];
    wki = c[j1];
    xr = a[j2 + 0] - a[k2 + 0];
    xi = a[j2 + 1] + a[k2 + 1];
    yr = wkr * xr + wki * xi;
    yi = wkr * xi - wki * xr;
    a[j2 + 0] = a[j2 + 0] - yr;
    a[j2 + 1] = yi - a[j2 + 1];
    a[k2 + 0] = yr + a[k2 + 0];
    a[k2 + 1] = yi - a[k2 + 1];
  }
  a[65] = -a[65];
}
#endif

#if !defined(WEBRTC_ARCH_X86_FAMILY) && !defined(WEBRTC_HAS_NEON)
// Fixed-point integer acceleration for 128-point Ooura FFT.
// Eliminates ~595 soft-float libcalls per FFT on architectures without scalar hardware FP (such as Xtensa HiFi).

static const int32_t rdft_w_q30[64] = {
    1073741824, 0, 759250112, 759250112,
    992008128, 410903232, 410903232, 992008128,
    1053110144, 209476640, 596539008, 892783680,
    892783680, 596539008, 209476640, 1053110144,
    1068571456, 105245104, 681174656, 830013632,
    946955776, 506158400, 311690784, 1027506880,
    1027506880, 311690784, 506158400, 946955776,
    830013632, 681174656, 105245104, 1068571456,
    759250112, 536224224, 534285728, 531060096,
    526555072, 520781568, 513753440, 505487616,
    496004064, 485325536, 473477888, 460489536,
    446391840, 431218752, 415006816, 397795136,
    379625056, 360540480, 340587328, 319813632,
    298269504, 276006816, 253079200, 229541888,
    205451616, 180866368, 155845392, 130449000,
    104738320, 78775320, 52622552, 26343008,
};

static const int32_t rdft_wk3ri_first_q30[16] = {
    1073741824, 0, 410903232, 992008128,
    892783616, 596539008, -209476671, 1053110144,
    1027506880, 311690816, 105245120, 1068571456,
    681174656, 830013696, -506158528, 946955648,
};

static const int32_t rdft_wk3ri_second_q30[16] = {
    -759250112, 759250112, -992008128, -410903232,
    -1053110144, 209476671, -596539008, -892783616,
    -946955648, 506158528, -830013696, -681174656,
    -1068571456, -105245120, -311690816, -1027506880,
};

static inline int32_t q30_mul(int32_t a, int32_t w_q30) {
  return (int32_t)(((int64_t)a * w_q30) >> 30);
}

static void bitrv2_128_int(int32_t* a) {
  unsigned int j, j1, k, k1;
  int32_t xr, xi, yr, yi;
  const int ip[4] = {0, 64, 32, 96};
  for (k = 0; k < 4; k++) {
    for (j = 0; j < k; j++) {
      j1 = 2 * j + ip[k];
      k1 = 2 * k + ip[j];
      xr = a[j1 + 0]; xi = a[j1 + 1];
      yr = a[k1 + 0]; yi = a[k1 + 1];
      a[j1 + 0] = yr; a[j1 + 1] = yi;
      a[k1 + 0] = xr; a[k1 + 1] = xi;
      j1 += 8; k1 += 16;
      xr = a[j1 + 0]; xi = a[j1 + 1];
      yr = a[k1 + 0]; yi = a[k1 + 1];
      a[j1 + 0] = yr; a[j1 + 1] = yi;
      a[k1 + 0] = xr; a[k1 + 1] = xi;
      j1 += 8; k1 -= 8;
      xr = a[j1 + 0]; xi = a[j1 + 1];
      yr = a[k1 + 0]; yi = a[k1 + 1];
      a[j1 + 0] = yr; a[j1 + 1] = yi;
      a[k1 + 0] = xr; a[k1 + 1] = xi;
      j1 += 8; k1 += 16;
      xr = a[j1 + 0]; xi = a[j1 + 1];
      yr = a[k1 + 0]; yi = a[k1 + 1];
      a[j1 + 0] = yr; a[j1 + 1] = yi;
      a[k1 + 0] = xr; a[k1 + 1] = xi;
    }
    j1 = 2 * k + 8 + ip[k];
    k1 = j1 + 8;
    xr = a[j1 + 0]; xi = a[j1 + 1];
    yr = a[k1 + 0]; yi = a[k1 + 1];
    a[j1 + 0] = yr; a[j1 + 1] = yi;
    a[k1 + 0] = xr; a[k1 + 1] = xi;
  }
}

static void cft1st_128_int(int32_t* a) {
  const int n = 128;
  int j, k1, k2;
  int32_t wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
  int32_t x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

  x0r = a[0] + a[2]; x0i = a[1] + a[3];
  x1r = a[0] - a[2]; x1i = a[1] - a[3];
  x2r = a[4] + a[6]; x2i = a[5] + a[7];
  x3r = a[4] - a[6]; x3i = a[5] - a[7];
  a[0] = x0r + x2r; a[1] = x0i + x2i;
  a[4] = x0r - x2r; a[5] = x0i - x2i;
  a[2] = x1r - x3i; a[3] = x1i + x3r;
  a[6] = x1r + x3i; a[7] = x1i - x3r;

  wk1r = rdft_w_q30[2];
  x0r = a[8] + a[10]; x0i = a[9] + a[11];
  x1r = a[8] - a[10]; x1i = a[9] - a[11];
  x2r = a[12] + a[14]; x2i = a[13] + a[15];
  x3r = a[12] - a[14]; x3i = a[13] - a[15];
  a[8] = x0r + x2r; a[9] = x0i + x2i;
  a[12] = x2i - x0i; a[13] = x0r - x2r;
  x0r = x1r - x3i; x0i = x1i + x3r;
  a[10] = q30_mul(x0r - x0i, wk1r);
  a[11] = q30_mul(x0r + x0i, wk1r);
  x0r = x3i + x1r; x0i = x3r - x1i;
  a[14] = q30_mul(x0i - x0r, wk1r);
  a[15] = q30_mul(x0i + x0r, wk1r);

  k1 = 0;
  for (j = 16; j < n; j += 16) {
    k1 += 2; k2 = 2 * k1;
    wk2r = rdft_w_q30[k1 + 0]; wk2i = rdft_w_q30[k1 + 1];
    wk1r = rdft_w_q30[k2 + 0]; wk1i = rdft_w_q30[k2 + 1];
    wk3r = rdft_wk3ri_first_q30[k1 + 0]; wk3i = rdft_wk3ri_first_q30[k1 + 1];
    x0r = a[j + 0] + a[j + 2]; x0i = a[j + 1] + a[j + 3];
    x1r = a[j + 0] - a[j + 2]; x1i = a[j + 1] - a[j + 3];
    x2r = a[j + 4] + a[j + 6]; x2i = a[j + 5] + a[j + 7];
    x3r = a[j + 4] - a[j + 6]; x3i = a[j + 5] - a[j + 7];
    a[j + 0] = x0r + x2r; a[j + 1] = x0i + x2i;
    x0r -= x2r; x0i -= x2i;
    a[j + 4] = q30_mul(x0r, wk2r) - q30_mul(x0i, wk2i);
    a[j + 5] = q30_mul(x0i, wk2r) + q30_mul(x0r, wk2i);
    x0r = x1r - x3i; x0i = x1i + x3r;
    a[j + 2] = q30_mul(x0r, wk1r) - q30_mul(x0i, wk1i);
    a[j + 3] = q30_mul(x0i, wk1r) + q30_mul(x0r, wk1i);
    x0r = x1r + x3i; x0i = x1i - x3r;
    a[j + 6] = q30_mul(x0r, wk3r) - q30_mul(x0i, wk3i);
    a[j + 7] = q30_mul(x0i, wk3r) + q30_mul(x0r, wk3i);
    wk1r = rdft_w_q30[k2 + 2]; wk1i = rdft_w_q30[k2 + 3];
    wk3r = rdft_wk3ri_second_q30[k1 + 0]; wk3i = rdft_wk3ri_second_q30[k1 + 1];
    x0r = a[j + 8] + a[j + 10]; x0i = a[j + 9] + a[j + 11];
    x1r = a[j + 8] - a[j + 10]; x1i = a[j + 9] - a[j + 11];
    x2r = a[j + 12] + a[j + 14]; x2i = a[j + 13] + a[j + 15];
    x3r = a[j + 12] - a[j + 14]; x3i = a[j + 13] - a[j + 15];
    a[j + 8] = x0r + x2r; a[j + 9] = x0i + x2i;
    x0r -= x2r; x0i -= x2i;
    a[j + 12] = -q30_mul(x0r, wk2i) - q30_mul(x0i, wk2r);
    a[j + 13] = -q30_mul(x0i, wk2i) + q30_mul(x0r, wk2r);
    x0r = x1r - x3i; x0i = x1i + x3r;
    a[j + 10] = q30_mul(x0r, wk1r) - q30_mul(x0i, wk1i);
    a[j + 11] = q30_mul(x0i, wk1r) + q30_mul(x0r, wk1i);
    x0r = x1r + x3i; x0i = x1i - x3r;
    a[j + 14] = q30_mul(x0r, wk3r) - q30_mul(x0i, wk3i);
    a[j + 15] = q30_mul(x0i, wk3r) + q30_mul(x0r, wk3i);
  }
}

static void cftmdl_128_int(int32_t* a) {
  const int l = 8;
  const int n = 128;
  const int m = 32;
  int j0, j1, j2, j3, k, k1, k2, m2;
  int32_t wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
  int32_t x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

  for (j0 = 0; j0 < l; j0 += 2) {
    j1 = j0 + 8; j2 = j0 + 16; j3 = j0 + 24;
    x0r = a[j0 + 0] + a[j1 + 0]; x0i = a[j0 + 1] + a[j1 + 1];
    x1r = a[j0 + 0] - a[j1 + 0]; x1i = a[j0 + 1] - a[j1 + 1];
    x2r = a[j2 + 0] + a[j3 + 0]; x2i = a[j2 + 1] + a[j3 + 1];
    x3r = a[j2 + 0] - a[j3 + 0]; x3i = a[j2 + 1] - a[j3 + 1];
    a[j0 + 0] = x0r + x2r; a[j0 + 1] = x0i + x2i;
    a[j2 + 0] = x0r - x2r; a[j2 + 1] = x0i - x2i;
    a[j1 + 0] = x1r - x3i; a[j1 + 1] = x1i + x3r;
    a[j3 + 0] = x1r + x3i; a[j3 + 1] = x1i - x3r;
  }
  wk1r = rdft_w_q30[2];
  for (j0 = m; j0 < l + m; j0 += 2) {
    j1 = j0 + 8; j2 = j0 + 16; j3 = j0 + 24;
    x0r = a[j0 + 0] + a[j1 + 0]; x0i = a[j0 + 1] + a[j1 + 1];
    x1r = a[j0 + 0] - a[j1 + 0]; x1i = a[j0 + 1] - a[j1 + 1];
    x2r = a[j2 + 0] + a[j3 + 0]; x2i = a[j2 + 1] + a[j3 + 1];
    x3r = a[j2 + 0] - a[j3 + 0]; x3i = a[j2 + 1] - a[j3 + 1];
    a[j0 + 0] = x0r + x2r; a[j0 + 1] = x0i + x2i;
    a[j2 + 0] = x2i - x0i; a[j2 + 1] = x0r - x2r;
    x0r = x1r - x3i; x0i = x1i + x3r;
    a[j1 + 0] = q30_mul(x0r - x0i, wk1r);
    a[j1 + 1] = q30_mul(x0r + x0i, wk1r);
    x0r = x3i + x1r; x0i = x3r - x1i;
    a[j3 + 0] = q30_mul(x0i - x0r, wk1r);
    a[j3 + 1] = q30_mul(x0i + x0r, wk1r);
  }
  k1 = 0; m2 = 2 * m;
  for (k = m2; k < n; k += m2) {
    k1 += 2; k2 = 2 * k1;
    wk2r = rdft_w_q30[k1 + 0]; wk2i = rdft_w_q30[k1 + 1];
    wk1r = rdft_w_q30[k2 + 0]; wk1i = rdft_w_q30[k2 + 1];
    wk3r = rdft_wk3ri_first_q30[k1 + 0]; wk3i = rdft_wk3ri_first_q30[k1 + 1];
    for (j0 = k; j0 < l + k; j0 += 2) {
      j1 = j0 + 8; j2 = j0 + 16; j3 = j0 + 24;
      x0r = a[j0 + 0] + a[j1 + 0]; x0i = a[j0 + 1] + a[j1 + 1];
      x1r = a[j0 + 0] - a[j1 + 0]; x1i = a[j0 + 1] - a[j1 + 1];
      x2r = a[j2 + 0] + a[j3 + 0]; x2i = a[j2 + 1] + a[j3 + 1];
      x3r = a[j2 + 0] - a[j3 + 0]; x3i = a[j2 + 1] - a[j3 + 1];
      a[j0 + 0] = x0r + x2r; a[j0 + 1] = x0i + x2i;
      x0r -= x2r; x0i -= x2i;
      a[j2 + 0] = q30_mul(x0r, wk2r) - q30_mul(x0i, wk2i);
      a[j2 + 1] = q30_mul(x0i, wk2r) + q30_mul(x0r, wk2i);
      x0r = x1r - x3i; x0i = x1i + x3r;
      a[j1 + 0] = q30_mul(x0r, wk1r) - q30_mul(x0i, wk1i);
      a[j1 + 1] = q30_mul(x0i, wk1r) + q30_mul(x0r, wk1i);
      x0r = x1r + x3i; x0i = x1i - x3r;
      a[j3 + 0] = q30_mul(x0r, wk3r) - q30_mul(x0i, wk3i);
      a[j3 + 1] = q30_mul(x0i, wk3r) + q30_mul(x0r, wk3i);
    }
    wk1r = rdft_w_q30[k2 + 2]; wk1i = rdft_w_q30[k2 + 3];
    wk3r = rdft_wk3ri_second_q30[k1 + 0]; wk3i = rdft_wk3ri_second_q30[k1 + 1];
    for (j0 = k + m; j0 < l + (k + m); j0 += 2) {
      j1 = j0 + 8; j2 = j0 + 16; j3 = j0 + 24;
      x0r = a[j0 + 0] + a[j1 + 0]; x0i = a[j0 + 1] + a[j1 + 1];
      x1r = a[j0 + 0] - a[j1 + 0]; x1i = a[j0 + 1] - a[j1 + 1];
      x2r = a[j2 + 0] + a[j3 + 0]; x2i = a[j2 + 1] + a[j3 + 1];
      x3r = a[j2 + 0] - a[j3 + 0]; x3i = a[j2 + 1] - a[j3 + 1];
      a[j0 + 0] = x0r + x2r; a[j0 + 1] = x0i + x2i;
      x0r -= x2r; x0i -= x2i;
      a[j2 + 0] = -q30_mul(x0r, wk2i) - q30_mul(x0i, wk2r);
      a[j2 + 1] = -q30_mul(x0i, wk2i) + q30_mul(x0r, wk2r);
      x0r = x1r - x3i; x0i = x1i + x3r;
      a[j1 + 0] = q30_mul(x0r, wk1r) - q30_mul(x0i, wk1i);
      a[j1 + 1] = q30_mul(x0i, wk1r) + q30_mul(x0r, wk1i);
      x0r = x1r + x3i; x0i = x1i - x3r;
      a[j3 + 0] = q30_mul(x0r, wk3r) - q30_mul(x0i, wk3i);
      a[j3 + 1] = q30_mul(x0i, wk3r) + q30_mul(x0r, wk3i);
    }
  }
}

static void cftfsub_128_int(int32_t* a) {
  int j, j1, j2, j3, l;
  int32_t x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

  cft1st_128_int(a);
  cftmdl_128_int(a);
  l = 32;
  for (j = 0; j < l; j += 2) {
    j1 = j + l; j2 = j1 + l; j3 = j2 + l;
    x0r = a[j] + a[j1]; x0i = a[j + 1] + a[j1 + 1];
    x1r = a[j] - a[j1]; x1i = a[j + 1] - a[j1 + 1];
    x2r = a[j2 + 0] + a[j3 + 0]; x2i = a[j2 + 1] + a[j3 + 1];
    x3r = a[j2 + 0] - a[j3 + 0]; x3i = a[j2 + 1] - a[j3 + 1];
    a[j] = x0r + x2r; a[j + 1] = x0i + x2i;
    a[j2] = x0r - x2r; a[j2 + 1] = x0i - x2i;
    a[j1] = x1r - x3i; a[j1 + 1] = x1i + x3r;
    a[j3] = x1r + x3i; a[j3 + 1] = x1i - x3r;
  }
}

static void cftbsub_128_int(int32_t* a) {
  int j, j1, j2, j3, l;
  int32_t x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

  cft1st_128_int(a);
  cftmdl_128_int(a);
  l = 32;
  for (j = 0; j < l; j += 2) {
    j1 = j + l; j2 = j1 + l; j3 = j2 + l;
    x0r = a[j] + a[j1];
    x0i = -a[j + 1] - a[j1 + 1];
    x1r = a[j] - a[j1];
    x1i = -a[j + 1] + a[j1 + 1];
    x2r = a[j2] + a[j3];
    x2i = a[j2 + 1] + a[j3 + 1];
    x3r = a[j2] - a[j3];
    x3i = a[j2 + 1] - a[j3 + 1];
    a[j] = x0r + x2r;
    a[j + 1] = x0i - x2i;
    a[j2] = x0r - x2r;
    a[j2 + 1] = x0i + x2i;
    a[j1] = x1r - x3i;
    a[j1 + 1] = x1i - x3r;
    a[j3] = x1r + x3i;
    a[j3 + 1] = x1i + x3r;
  }
}

static void rftfsub_128_int(int32_t* a) {
  const int32_t* c = rdft_w_q30 + 32;
  int j1, j2, k1, k2;
  int32_t wkr, wki, xr, xi, yr, yi;

  for (j1 = 1, j2 = 2; j2 < 64; j1 += 1, j2 += 2) {
    k2 = 128 - j2;
    k1 = 32 - j1;
    wkr = (1 << 29) - c[k1];
    wki = c[j1];
    xr = a[j2 + 0] - a[k2 + 0];
    xi = a[j2 + 1] + a[k2 + 1];
    yr = q30_mul(wkr, xr) - q30_mul(wki, xi);
    yi = q30_mul(wkr, xi) + q30_mul(wki, xr);
    a[j2 + 0] -= yr;
    a[j2 + 1] -= yi;
    a[k2 + 0] += yr;
    a[k2 + 1] -= yi;
  }
}

static void rftbsub_128_int(int32_t* a) {
  const int32_t* c = rdft_w_q30 + 32;
  int j1, j2, k1, k2;
  int32_t wkr, wki, xr, xi, yr, yi;

  a[1] = -a[1];
  for (j1 = 1, j2 = 2; j2 < 64; j1 += 1, j2 += 2) {
    k2 = 128 - j2;
    k1 = 32 - j1;
    wkr = (1 << 29) - c[k1];
    wki = c[j1];
    xr = a[j2 + 0] - a[k2 + 0];
    xi = a[j2 + 1] + a[k2 + 1];
    yr = q30_mul(wkr, xr) + q30_mul(wki, xi);
    yi = q30_mul(wkr, xi) - q30_mul(wki, xr);
    a[j2 + 0] = a[j2 + 0] - yr;
    a[j2 + 1] = yi - a[j2 + 1];
    a[k2 + 0] = yr + a[k2 + 0];
    a[k2 + 1] = yi - a[k2 + 1];
  }
  a[65] = -a[65];
}

static inline int32_t float_bits_to_scaled_int(uint32_t u, int max_exp) {
  int exp = (u >> 23) & 0xff;
  if (exp == 0) return 0;

  int diff = max_exp - exp;
  if (diff > 25) return 0;

  int32_t mantissa = (1 << 23) | (u & 0x7fffff);
  int32_t val = (diff == 0) ? (mantissa << 1) : (mantissa >> (diff - 1));
  return ((int32_t)u < 0) ? -val : val;
}

static inline uint32_t scaled_int_to_float_bits(int32_t val, int max_exp) {
  if (val == 0) return 0;
  uint32_t sign = 0;
  if (val < 0) {
    sign = 0x80000000;
    val = -val;
  }
  int lz = __builtin_clz(val);
  int shift = 31 - lz - 23;
  uint32_t mantissa = (shift >= 0) ? (((uint32_t)val >> shift) & 0x7fffff)
                                   : (((uint32_t)val << (-shift)) & 0x7fffff);
  int exp = max_exp + 7 - lz;
  if (exp <= 0) return 0;
  if (exp >= 255) exp = 254;

  return sign | ((uint32_t)exp << 23) | mantissa;
}

static void FastFftInt(float* a) {
  uint32_t* u_a = reinterpret_cast<uint32_t*>(a);
  uint32_t max_u = 0;
  for (int i = 0; i < 128; i += 4) {
    uint32_t u0 = u_a[i] & 0x7fffffff;
    uint32_t u1 = u_a[i + 1] & 0x7fffffff;
    uint32_t u2 = u_a[i + 2] & 0x7fffffff;
    uint32_t u3 = u_a[i + 3] & 0x7fffffff;
    if (u0 > max_u) max_u = u0;
    if (u1 > max_u) max_u = u1;
    if (u2 > max_u) max_u = u2;
    if (u3 > max_u) max_u = u3;
  }
  if (max_u < 0x2b800000) {
    memset(a, 0, 128 * sizeof(float));
    return;
  }
  int max_exp = (max_u >> 23) & 0xff;

  int32_t a_q[128];
  for (int i = 0; i < 128; ++i) {
    a_q[i] = float_bits_to_scaled_int(u_a[i], max_exp);
  }

  bitrv2_128_int(a_q);
  cftfsub_128_int(a_q);
  rftfsub_128_int(a_q);
  int32_t xi = a_q[0] - a_q[1];
  a_q[0] += a_q[1];
  a_q[1] = xi;

  for (int i = 0; i < 128; ++i) {
    u_a[i] = scaled_int_to_float_bits(a_q[i], max_exp);
  }
}

static void FastInverseFftInt(float* a) {
  uint32_t* u_a = reinterpret_cast<uint32_t*>(a);
  uint32_t max_u = 0;
  for (int i = 0; i < 128; i += 4) {
    uint32_t u0 = u_a[i] & 0x7fffffff;
    uint32_t u1 = u_a[i + 1] & 0x7fffffff;
    uint32_t u2 = u_a[i + 2] & 0x7fffffff;
    uint32_t u3 = u_a[i + 3] & 0x7fffffff;
    if (u0 > max_u) max_u = u0;
    if (u1 > max_u) max_u = u1;
    if (u2 > max_u) max_u = u2;
    if (u3 > max_u) max_u = u3;
  }
  if (max_u < 0x2b800000) {
    memset(a, 0, 128 * sizeof(float));
    return;
  }
  int max_exp = (max_u >> 23) & 0xff;

  int32_t a_q[128];
  for (int i = 0; i < 128; ++i) {
    a_q[i] = float_bits_to_scaled_int(u_a[i], max_exp);
  }

  int32_t a1 = (a_q[0] - a_q[1]) >> 1;
  a_q[0] -= a1;
  a_q[1] = a1;

  rftbsub_128_int(a_q);
  bitrv2_128_int(a_q);
  cftbsub_128_int(a_q);

  for (int i = 0; i < 128; ++i) {
    u_a[i] = scaled_int_to_float_bits(a_q[i], max_exp);
  }
}
#endif

}  // namespace

OouraFft::OouraFft(bool sse2_available) {
#if defined(WEBRTC_ARCH_X86_FAMILY)
  use_sse2_ = sse2_available;
#else
  use_sse2_ = false;
#endif
}

OouraFft::OouraFft() {
#if defined(WEBRTC_ARCH_X86_FAMILY)
  use_sse2_ = (GetCPUInfo(kSSE2) != 0);
#else
  use_sse2_ = false;
#endif
}

OouraFft::~OouraFft() = default;

void OouraFft::Fft(float* a) const {
#if !defined(WEBRTC_ARCH_X86_FAMILY) && !defined(WEBRTC_HAS_NEON)
  FastFftInt(a);
#else
  float xi;
  bitrv2_128(a);
  cftfsub_128(a);
  rftfsub_128(a);
  xi = a[0] - a[1];
  a[0] += a[1];
  a[1] = xi;
#endif
}
void OouraFft::InverseFft(float* a) const {
#if !defined(WEBRTC_ARCH_X86_FAMILY) && !defined(WEBRTC_HAS_NEON)
  FastInverseFftInt(a);
#else
  a[1] = 0.5f * (a[0] - a[1]);
  a[0] -= a[1];
  rftbsub_128(a);
  bitrv2_128(a);
  cftbsub_128(a);
#endif
}

void OouraFft::cft1st_128(float* a) const {

#if defined(MIPS_FPU_LE)
  cft1st_128_mips(a);
#elif defined(WEBRTC_HAS_NEON)
  cft1st_128_neon(a);
#elif defined(WEBRTC_ARCH_X86_FAMILY)
  if (use_sse2_) {
    cft1st_128_SSE2(a);
  } else {
    cft1st_128_C(a);
  }
#else
  cft1st_128_C(a);
#endif
}
void OouraFft::cftmdl_128(float* a) const {
#if defined(MIPS_FPU_LE)
  cftmdl_128_mips(a);
#elif defined(WEBRTC_HAS_NEON)
  cftmdl_128_neon(a);
#elif defined(WEBRTC_ARCH_X86_FAMILY)
  if (use_sse2_) {
    cftmdl_128_SSE2(a);
  } else {
    cftmdl_128_C(a);
  }
#else
  cftmdl_128_C(a);
#endif
}
void OouraFft::rftfsub_128(float* a) const {
#if defined(MIPS_FPU_LE)
  rftfsub_128_mips(a);
#elif defined(WEBRTC_HAS_NEON)
  rftfsub_128_neon(a);
#elif defined(WEBRTC_ARCH_X86_FAMILY)
  if (use_sse2_) {
    rftfsub_128_SSE2(a);
  } else {
    rftfsub_128_C(a);
  }
#else
  rftfsub_128_C(a);
#endif
}

void OouraFft::rftbsub_128(float* a) const {
#if defined(MIPS_FPU_LE)
  rftbsub_128_mips(a);
#elif defined(WEBRTC_HAS_NEON)
  rftbsub_128_neon(a);
#elif defined(WEBRTC_ARCH_X86_FAMILY)
  if (use_sse2_) {
    rftbsub_128_SSE2(a);
  } else {
    rftbsub_128_C(a);
  }
#else
  rftbsub_128_C(a);
#endif
}

void OouraFft::cftbsub_128(float* a) const {
  int j, j1, j2, j3, l;
  float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

  cft1st_128(a);
  cftmdl_128(a);
  l = 32;

  for (j = 0; j < l; j += 2) {
    j1 = j + l;
    j2 = j1 + l;
    j3 = j2 + l;
    x0r = a[j] + a[j1];
    x0i = -a[j + 1] - a[j1 + 1];
    x1r = a[j] - a[j1];
    x1i = -a[j + 1] + a[j1 + 1];
    x2r = a[j2] + a[j3];
    x2i = a[j2 + 1] + a[j3 + 1];
    x3r = a[j2] - a[j3];
    x3i = a[j2 + 1] - a[j3 + 1];
    a[j] = x0r + x2r;
    a[j + 1] = x0i - x2i;
    a[j2] = x0r - x2r;
    a[j2 + 1] = x0i + x2i;
    a[j1] = x1r - x3i;
    a[j1 + 1] = x1i - x3r;
    a[j3] = x1r + x3i;
    a[j3 + 1] = x1i + x3r;
  }
}

void OouraFft::cftfsub_128(float* a) const {
  int j, j1, j2, j3, l;
  float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

  cft1st_128(a);
  cftmdl_128(a);
  l = 32;
  for (j = 0; j < l; j += 2) {
    j1 = j + l;
    j2 = j1 + l;
    j3 = j2 + l;
    x0r = a[j] + a[j1];
    x0i = a[j + 1] + a[j1 + 1];
    x1r = a[j] - a[j1];
    x1i = a[j + 1] - a[j1 + 1];
    x2r = a[j2] + a[j3];
    x2i = a[j2 + 1] + a[j3 + 1];
    x3r = a[j2] - a[j3];
    x3i = a[j2 + 1] - a[j3 + 1];
    a[j] = x0r + x2r;
    a[j + 1] = x0i + x2i;
    a[j2] = x0r - x2r;
    a[j2 + 1] = x0i - x2i;
    a[j1] = x1r - x3i;
    a[j1 + 1] = x1i + x3r;
    a[j3] = x1r + x3i;
    a[j3 + 1] = x1i - x3r;
  }
}

void OouraFft::bitrv2_128(float* a) const {
  /*
      Following things have been attempted but are no faster:
      (a) Storing the swap indexes in a LUT (index calculations are done
          for 'free' while waiting on memory/L1).
      (b) Consolidate the load/store of two consecutive floats by a 64 bit
          integer (execution is memory/L1 bound).
      (c) Do a mix of floats and 64 bit integer to maximize register
          utilization (execution is memory/L1 bound).
      (d) Replacing ip[i] by ((k<<31)>>25) + ((k >> 1)<<5).
      (e) Hard-coding of the offsets to completely eliminates index
          calculations.
  */

  unsigned int j, j1, k, k1;
  float xr, xi, yr, yi;

  const int ip[4] = {0, 64, 32, 96};
  for (k = 0; k < 4; k++) {
    for (j = 0; j < k; j++) {
      j1 = 2 * j + ip[k];
      k1 = 2 * k + ip[j];
      xr = a[j1 + 0];
      xi = a[j1 + 1];
      yr = a[k1 + 0];
      yi = a[k1 + 1];
      a[j1 + 0] = yr;
      a[j1 + 1] = yi;
      a[k1 + 0] = xr;
      a[k1 + 1] = xi;
      j1 += 8;
      k1 += 16;
      xr = a[j1 + 0];
      xi = a[j1 + 1];
      yr = a[k1 + 0];
      yi = a[k1 + 1];
      a[j1 + 0] = yr;
      a[j1 + 1] = yi;
      a[k1 + 0] = xr;
      a[k1 + 1] = xi;
      j1 += 8;
      k1 -= 8;
      xr = a[j1 + 0];
      xi = a[j1 + 1];
      yr = a[k1 + 0];
      yi = a[k1 + 1];
      a[j1 + 0] = yr;
      a[j1 + 1] = yi;
      a[k1 + 0] = xr;
      a[k1 + 1] = xi;
      j1 += 8;
      k1 += 16;
      xr = a[j1 + 0];
      xi = a[j1 + 1];
      yr = a[k1 + 0];
      yi = a[k1 + 1];
      a[j1 + 0] = yr;
      a[j1 + 1] = yi;
      a[k1 + 0] = xr;
      a[k1 + 1] = xi;
    }
    j1 = 2 * k + 8 + ip[k];
    k1 = j1 + 8;
    xr = a[j1 + 0];
    xi = a[j1 + 1];
    yr = a[k1 + 0];
    yi = a[k1 + 1];
    a[j1 + 0] = yr;
    a[j1 + 1] = yi;
    a[k1 + 0] = xr;
    a[k1 + 1] = xi;
  }
}

}  // namespace webrtc
