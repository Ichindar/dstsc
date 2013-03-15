/*
------------------------------------------------------------------
 Description: Generates an (Lempel-Ziv) factorisation using the 
              Ohlebusch-Cog LZ algorithm (lzOG). Algorithm uses
	      suffix array to LZ factorise string.
	      C++ implementation (cpp) file.
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

#include "lzOG.h"

void sop(LONGINT i, LONGINT l, LONGINT j, LONGINT *lps, LONGINT *prev_occurance, int null_value)
    { 
   if(j == 0 and l ==0 and i==0)
       return;
   // assert(i>j); Used for debugging
   // assert(i!=j);
   if(lps[i] == null_value)
       {
       lps[i] = l;
       prev_occurance[i] = j;
       }
   else
       {
       if(lps[i] < l)
	   {
           if(prev_occurance[i] > j)  // if the ith previous occurance is smaller than j, swap them 
	       sop(prev_occurance[i], lps[i], j, lps, prev_occurance, null_value);
	   else
	       sop(j, lps[i],prev_occurance[i], lps, prev_occurance, null_value);
	   lps[i] = l;
	   prev_occurance[i] = j;
	   }
       else
	   {
           if(prev_occurance[i] > j)  // if the ith previous occurance is smaller than j, swap them  
	       sop(prev_occurance[i], l, j, lps, prev_occurance, null_value);
	   else
	       sop(j, l, prev_occurance[i],lps, prev_occurance, null_value);			
	   }
	}
}

void lz_factorise_file(const char * filename, vector<LONGINT> &offsets, vector<LONGINT> &lengths)
{
/* 
Procedure loads and factorises 'filename' file. 
Returns two vectors one with the offsets and one with the lengths
*/

    GtUchar *filecontents; // File contents will be stored here
    LONGINT numofbytes;
    filecontents = (GtUchar*)creatememorymap(filename, &numofbytes);
    lz_factorise(numofbytes, filecontents, offsets, lengths);
    deletememorymap(filecontents, numofbytes);
}

void lz_factorise(LONGINT bytenum, const unsigned char * filecontents, vector<LONGINT> &offsets, vector<LONGINT> &lengths)
{ // procedure accepts suffix array of byte array 'contents' of size 'bytenum'

  // Step 1: Compute Suffix Array

   LONGINT *sa;           // Suffix array will be stored here
   sa = gt_sain_plain_sortsuffixes(filecontents, bytenum, false);

  // Step 2: Compute phi array

  LONGINT *phi = new LONGINT[bytenum];
  LONGINT *previous_occurance = new LONGINT[bytenum];
  
  for(int i=0; i<bytenum; i++)
     previous_occurance[i] = filecontents[i];

  phi[sa[0]] = sa[bytenum-1];
  for(int i=1; i<bytenum; i++)
     phi[sa[i]] = sa[i-1];

  // Step 3: Calculate Longest Previous Substrings (LPS) and Previous Occurances

  for(int i=0; i<bytenum; ++i)
      sa[i] = -1;  // To save space, we reuse the SA array to store LPS values
  
  int maximumFactorLength = 0;
  int l = 0;
  int j;
  for(int i=0; i<bytenum; i++)
      {  
      j = phi[i];  
      while(filecontents[i+l] == filecontents[j+l])
	  l++;
      if (i>j)
	  sop(i, l, j, sa, previous_occurance, -1);
      else
	  sop(j, l, i, sa, previous_occurance, -1);
      if (l > 0)
	  l--;
      }

  // Step 4: Perform LZ factorisation using LPS and Previous Occurances

  int i = 1;

  offsets.push_back(filecontents[0]); lengths.push_back(0); // adding first byte
  
  while(i < bytenum)
      {
      if(sa[i] == 0)
	  {
	  offsets.push_back(filecontents[i]); 
	  lengths.push_back(sa[i]);
	  }
      else
	  {
	  offsets.push_back(previous_occurance[i]);
	  lengths.push_back(sa[i]);
	  }
      if (sa[i] > 1)
	  i = i + sa[i];
      else
	  i++;
      }

  // Finalisation: Deallocate memory
  gt_free(sa);
  delete[] phi;
  delete[] previous_occurance;
}
