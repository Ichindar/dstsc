/*
------------------------------------------------------------------
 Description: Simple program illustrating the use of the lzOG
              Lempel_Ziv factorisation library.
 Author: Angelos Molfetas (2013)
 Copyright: The University of Melbourne (2013)
 Licence: BSD licence, see attached LICENCE file
 -----------------------------------------------------------------
*/

#include "lzOG.h"
#include <vector>
#include <iostream>

using namespace std;

int main(int argc, const char **argv)
{


   LONGINT numofbytes, *suftab;
   vector<LONGINT> offsets; // Store the factor offsets
   vector<LONGINT> lengths; // Store the length offsets

   if (argc != 2) // Validate input
       {
       fprintf(stderr,"Usage: %s <inputfile>\n",argv[0]);
       exit(EXIT_FAILURE);
       }

   lz_factorise_file(argv[1], offsets, lengths);
   
   for(int i = 0; i < offsets.size() - 1; i++) // Print the lz factors
       if (lengths[i] != 0)
	   cout << "(" << offsets[i] << "," << lengths[i] << ")" << endl;
       else  // when length is zero we are storing a char instead of an offset
	   cout << "(" << (unsigned char)offsets[i] << "," << lengths[i] << ")" << endl;

   exit(EXIT_SUCCESS);
}
