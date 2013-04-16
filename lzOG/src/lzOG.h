/*
------------------------------------------------------------------
 Description: Generates an (Lempel-Ziv) factorisation using suffix 
              array using the Ohlebusch-Cog LZ algorithm (lzOG).  
	      C++ header (h) file.
 Author:  Angelos Molfetas (2013)
 Credits: Thanks to Simon Cog. LZ factorisor method largerly based 
          on implementation he provided.
	  Thanks to Stefan Kurtz. Implementation relies Suffix
          Array implementation by Stefan Kurtz.
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
#include <vector>
#include "../../SAIS-SK/src/sk-sain.h"
#include "../../SAIS-SK/src/gt-alloc.h"
#include "../../SAIS-SK/src/mmap.h"

using namespace std;

typedef long unsigned int LONGINT;

void sop(LONGINT i, LONGINT l, LONGINT j, LONGINT *lps, LONGINT *prev_occurance, int null_value);
void lz_factorise(LONGINT bytenum, const unsigned char * filecontents, vector<LONGINT> &offsets, vector<LONGINT> &lengths);
void lz_refs(LONGINT bytenum, const unsigned char * filecontents, vector<LONGINT> &offsets, vector<LONGINT> &lengths);
void lz_factorise_file(const char* filename, vector<LONGINT> &offsets, vector<LONGINT> &lengths);
