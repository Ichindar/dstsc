/*
------------------------------------------------------------------
 Description: Simple program illustrating the use of the lzOG
              Lempel_Ziv factorisation library to factorise a
	      string.
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
   LONGINT *suftab;
   vector<LONGINT> offsets; // Store the factor offsets
   vector<LONGINT> lengths; // Store the length offsets

   const char* name = "123451234512345";

   lz_factorise(16, (unsigned char *) name, offsets, lengths);
   
   for(int i = 0; i < offsets.size(); i++) // Print the lz factors
       if (lengths[i] != 0)
	   cout << "(" << offsets[i] << "," << lengths[i] << ")" << endl;
       else  // when length is zero we are storing a char instead of an offset
	   cout << "(" << (unsigned char)offsets[i] << "," << lengths[i] << ")" << endl;

   exit(EXIT_SUCCESS);
}
