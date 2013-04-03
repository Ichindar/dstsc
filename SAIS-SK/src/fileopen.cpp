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
#include <sys/types.h>  // for fstat
#include <sys/stat.h>   // for fstat
#include <unistd.h>     // for fstat
#include <fcntl.h>      // for open
#include "fileopen.h"

int simplefileOpen(const char *filename,unsigned long *numofbytes)
{
  int fd;
  struct stat buf;

  fd = open(filename,O_RDONLY); /* open file for read, get filedes */
  if (fd == -1) /* check for error code */
  {
    fprintf(stderr,"Cannot open \"%s\"\n",filename);
    return -1;
  }
  if (fstat(fd,&buf) == -1)  /* get status information of file */
  {
    fprintf(stderr,"Cannot access status of file: %s\n",filename);
    return -2;
  }
  *numofbytes = (unsigned long)buf.st_size; /* store file size in address of numofbytes */
  return fd;
}
