/*
  Copyright (c) 2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

  Additional credits: 
  -- Angelos Molfetas (2013) Added explicit type casting for void 
     pointers so code compiles with C++ compilier.
*/

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include "sk-sain.h"
#include "gt-alloc.h"
#include "mmap.h"
#include "sufcheck.h"

int main(int argc,const char **argv)
{
  GtUchar *filecontents;
  unsigned long numofbytes, *suftab;

  if (argc != 2)
  {
    fprintf(stderr,"Usage: %s <inputfile>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  filecontents = (GtUchar *) creatememorymap(argv[1],&numofbytes);
  suftab = gt_sain_plain_sortsuffixes(filecontents,numofbytes,false);
  gt_suftab_lightweightcheck(filecontents,numofbytes,suftab);
  printf("[");
  for(int i = 0; i < numofbytes-1; i++)
  {
  	printf("%lu, ", suftab[i]);
  }
  printf("%lu]\n", suftab[numofbytes-1]);
  gt_free(suftab);
  deletememorymap(filecontents,numofbytes);
  exit(EXIT_SUCCESS);
}
