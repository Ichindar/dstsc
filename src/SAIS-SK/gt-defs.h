/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef GT_DEFS_H
#define GT_DEFS_H

#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>

typedef unsigned char GtUchar;

#ifdef _LP64

#define GT_LOGWORDSIZE    6         /* base 2 logarithm of wordsize */
typedef uint64_t GtBitsequence;
#else

#define GT_LOGWORDSIZE   5          /* base 2 logarithm of wordsize */
typedef uint32_t GtBitsequence;

#endif

#define GT_INTWORDSIZE\
        (1 << GT_LOGWORDSIZE) /* # of bits in unsigned long = w */

#define GT_FIRSTTWOBITS\
        (((GtBitsequence) 3) << (GT_INTWORDSIZE-2)) /* \(11^{w-2}\) */

#define GT_MULT2(N)     ((N) << 1)
#define GT_MULT4(N)     ((N) << 2)
#define GT_DIV2(N)      ((N) >> 1)

#define GT_UNUSED \
        __attribute__ ((unused)) /*@unused@*/

#define gt_assert(EXPR)  assert(EXPR)

#endif
