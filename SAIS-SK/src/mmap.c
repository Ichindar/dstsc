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

#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include "fileopen.h"

void *creatememorymap(const char *filename,
                      unsigned long *numofbytes)
{
  int fd;
  void *memorymap;

  fd = simplefileOpen(filename,numofbytes);
  if (fd < 0)
  {
    return NULL;
  }
  memorymap = mmap(NULL, /* no address specified for map */
                   (size_t) *numofbytes,
                   PROT_READ, /* pages may be read */
                   MAP_PRIVATE, /* unclear why need for reading */
                   fd, /* file descriptor */
                   (off_t) 0); /* offset: multiple of page size */
  if (memorymap == (void *) MAP_FAILED)
  {
    fprintf(stderr,"%s(%s) failed\n",__func__,filename);
    exit(EXIT_FAILURE);
  }
  return memorymap;
}

void deletememorymap(void *memorymap,unsigned long numofbytes)
{
  if (munmap(memorymap,(size_t) numofbytes) != 0)
  {
    fprintf(stderr,"%s failed\n",__func__);
    exit(EXIT_FAILURE);
  }
}
