/*
------------------------------------------------------------------
 Description: Generates an (Lempel-Ziv) factorisation using suffix 
              array. Uses algorithm by Enno Ohlebusch and Simon 
	      Cog (2011). Relies on Suffix Array implementation by
	      Stefan Kurtz.
 Author: Angelos Molfetas (2013)
 Copyright: The University of Melbourne (2013)
 Licence: BSD licence, see attached LICENCE file
 Reference: Ohlebusch, E., Cog, S. (2011). Combinatorial Pattern
            Matching 2011 (CPM 2011). Lecture Notes for Computer
	    Science. Volume 6661. pp. 15 - 26. Springer-Verlag 
	    Berlin.
 -----------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include "../SAIS-SK/sk-sain.h"
#include "../SAIS-SK/gt-alloc.h"
#include "../SAIS-SK/mmap.h"
//include "sufcheck.h"

int main(int argc, const char **argv)
{

   GtUchar *filecontents;
   unsigned long numofbytes, *suftab;
 
 
   // Validate input
   if (argc != 2)
   {
     fprintf(stderr,"Usage: %s <inputfile>\n",argv[0]);
     exit(EXIT_FAILURE);
   }

   filecontents = (GtUchar*)creatememorymap(argv[1],&numofbytes);
   suftab = gt_sain_plain_sortsuffixes(filecontents,numofbytes,false);
   gt_free(suftab);
   deletememorymap(filecontents,numofbytes);
   exit(EXIT_SUCCESS);
}
