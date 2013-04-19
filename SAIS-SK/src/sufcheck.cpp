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

  History: 
  -- Angelos Molfetas (2013) Added explicit type casting for void 
     pointers so code compiles with C++ compilier.
*/

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include "gt-alloc.h"
#include "gt-defs.h"

#define GT_DIVWORDSIZE(I)\
        ((I) >> GT_LOGWORDSIZE)              /* \((I) div w\) */

#define GT_FIRSTBIT\
        (((GtBitsequence) 1) << (GT_INTWORDSIZE-1)) /* \(10^{w-1}\) */

#define GT_ITHBIT(I)\
        (GT_FIRSTBIT >> (I))                 /* \(0^{i}10^{w-i-1}\)  */

/*
  \texttt{SETIBIT(TAB,I)} sets the \texttt{I}-th bit in bitarray
  \texttt{TAB} to 1.
*/

#define GT_SETIBIT(TAB,I)    (TAB)[GT_DIVWORDSIZE(I)] |= \
                                    GT_ITHBIT(GT_MODWORDSIZE(I))
/*
  \texttt{ISIBITSET(TAB,I)} checks if the \texttt{I}-th bit in bitarray
  \texttt{TAB} is 1.
*/

#define GT_ISIBITSET(TAB,I)  ((TAB)[GT_DIVWORDSIZE(I)] & \
                                    GT_ITHBIT(GT_MODWORDSIZE(I)))

#define GT_MODWORDSIZE(I)\
        ((I) & (GT_INTWORDSIZE-1))           /* \((I) mod w\) */

#define GT_NUMOFINTSFORBITS(N)\
        ((GT_DIVWORDSIZE(N) == 0)\
           ? (size_t) 1 \
           : ((size_t) 1 + (size_t) GT_DIVWORDSIZE((N) - 1)))

/*
  The following macro allocates a bitarray of \texttt{N} bits. All bits
  are off.
*/

#define GT_INITBITTAB(TAB,NUMOFBITS)\
        {\
          size_t tabsize = GT_NUMOFINTSFORBITS(NUMOFBITS);\
          TAB = (GtBitsequence*) gt_malloc(sizeof (GtBitsequence) * tabsize);\
          (void) memset(TAB,0,sizeof (GtBitsequence) * tabsize);\
        }

typedef struct
{
  unsigned long start, end;
  GtUchar firstchar;
} GtRangewithchar;

/* The following funktion implements the linear time algorithm of
   @INPROCEEDINGS{BUR:KAER:2003,
   author = {Burkhardt, S. and K{\"a}rkk{\"a}inen, J.},
   title = {{Fast Lightweight Suffix Array Construction and Checking}},
   booktitle = {{Proceedings of the 14th Annual Symposium on Combinatorial
                 Pattern Matching (CPM)}},
   year = {2003},
   editor = {{Baeza-Yates, R. and Ch{\'a}vez, E. and Crochemore, M.}},
   volume = {2676},
   series = {LNCS},
   pages = {200-210},
   publisher = {Springer-Verlag}
   }
   to check the following suffix-order condition of the sorted suffix array:
   For all characters c, if SA[i,j] contains the suffixes starting
   with charcter c, then SA[i]+1, SA[i+1]+1, \ldots, SA[j]+1 occur
   in SA in this order (but not consecutively in general).
   The running time of the algorithm is independent of the alphabet size.
   The main problem is that it requires random access to the sequence
   which slows it down.
*/

static void gt_suftab_bk_suffixorder(const GtUchar *sequence,
                                     unsigned long totallength,
                                     unsigned int numofchars,
                                     const unsigned long *suftab,
                                     const GtRangewithchar *rangestore,
                                     unsigned int numofranges)
{
  unsigned int rangeidx;
  unsigned long idx;
  unsigned long *nexttab = (unsigned long*) gt_calloc((size_t) numofchars,sizeof(*nexttab));

  for (rangeidx = 0; rangeidx < numofranges; rangeidx++)
  {
    nexttab[rangestore[rangeidx].firstchar] = rangestore[rangeidx].start;
  }
  for (idx = 0; idx < totallength; idx++)
  {
    unsigned long position = suftab[idx];

    if (position > 0)
    {
      GtUchar cc = sequence[position-1];
      unsigned long checkpos = suftab[nexttab[(int) cc]] + 1;
      if (checkpos != position)
      {
        fprintf(stderr,"idx=%lu,checkpos=%lu,position=%lu\n",
                        idx,checkpos,position);
        exit(EXIT_FAILURE);
      }
      nexttab[(int) cc]++;
    }
  }
  gt_free(nexttab);
}

void gt_suftab_lightweightcheck(const GtUchar *sequence,
                               unsigned long totallength,
                               const unsigned long *suftab)
{
  unsigned long idx, countbitsset = 0, previouspos = 0,
                firstspecial = totallength, rangestart = 0;
  const unsigned int numofchars = UCHAR_MAX;
  unsigned long charcount[UCHAR_MAX+1] = {0};
  unsigned int charidx, rangeidx = 0, numofranges;
  GtBitsequence *startposoccurs;
  GtUchar previouscc = 0;
  GtRangewithchar *rangestore;

  GT_INITBITTAB(startposoccurs,totallength+1);
  rangestore = (GtRangewithchar *) gt_malloc(sizeof(*rangestore) * numofchars);
  for (idx = 0; idx < totallength; idx++)
  {
    unsigned long position = suftab[idx];
    GtUchar cc;

    charcount[(int) sequence[idx]]++;
    if (GT_ISIBITSET(startposoccurs,position))
    {
      fprintf(stderr,"ERROR: suffix with startpos %lu already occurs\n",
              suftab[idx]);
      exit(EXIT_FAILURE);
    }
    GT_SETIBIT(startposoccurs,position);
    countbitsset++;
    cc = sequence[position];
    if (idx > 0)
    {
      if (previouscc > cc)
      {
        fprintf(stderr,"incorrect order: %lu = %lu=%u > %u=%lu=%lu\n",
                idx-1,position,(unsigned int) previouscc,
                (unsigned int) cc,previouspos,idx);
        exit(EXIT_FAILURE);
      } else
      {
        if (previouscc < cc)
        {
          gt_assert(rangeidx < numofchars);
          rangestore[rangeidx].start = rangestart;
          rangestore[rangeidx].end = idx-1;
          rangestore[rangeidx++].firstchar = previouscc;
          rangestart = idx;
        }
      }
    }
    previouscc = cc;
    previouspos = position;
  }
  if (countbitsset != totallength)
  {
    fprintf(stderr,"ERROR: only %lu of %lu suffixes occur\n",countbitsset,
                                totallength);
    exit(EXIT_FAILURE);
  }
  gt_free(startposoccurs);
  if (firstspecial == totallength)
  {
    gt_assert(firstspecial > 0 && rangeidx < numofchars);
    rangestore[rangeidx].start = rangestart;
    rangestore[rangeidx].end = firstspecial-1;
    rangestore[rangeidx++].firstchar = previouscc;
  }
  numofranges = rangeidx;
  for (charidx = 0, rangeidx = 0; charidx < numofchars; charidx++)
  {
    unsigned long count = charcount[charidx];
    if (count != 0)
    {
      gt_assert(rangestore[rangeidx].firstchar == (GtUchar) charidx);
      gt_assert(rangestore[rangeidx].end - rangestore[rangeidx].start + 1
                == count);
      rangeidx++;
    }
  }
  gt_suftab_bk_suffixorder(sequence,
                           totallength,
                           numofchars,
                           suftab,
                           rangestore,
                           numofranges);
  gt_free(rangestore);
}
