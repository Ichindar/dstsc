DSTSC - Data Structures and related algorithms
==============================================

C/C++ implementations for data structures and related algorithms.

Contents
--------

* Suffix Array construction algorithm 
* LZ factorisation algorithm

### Suffix Array construction

## Description

Takes a byte string and returns an array of integers corresponding to the sorted order of all suffixes in the byte string. 

For example, the suffixes of the the string "ABCAB" are as follows:

    Offset:   01234
    Suffix 1: ABCAB
    Suffix 2:  BCAB
    Suffix 3:   CAB
    Suffix 4:    AB
    Suffix 5:     B

The ordered suffixes for this example and their corresponding starting offsets are as follows:

    Offset:  012345
    Suffix 4:    AB  ---> starts at offset 3
    Suffix 1: ABCAB  ---> starts at offset 0
    Suffix 5:     B  ---> starts at offset 4
    Suffix 2:  BCAB  ---> starts at offset 1
    Suffix 3:   CAB  ---> starts at offset 2

The Suffix Array for string "ABCAB" is these starting offsets, that is:

    SA = [3, 0, 4, 1, 2]

## Implementation details

This implementation, SAIS_SK, is by Stefan Kurtz and uses the induction sorting algorithm (Nong, Zhang and Chan, 2011). Based on the ideas in Yuta Mori's implementation (http://code.google.com/p/libdivsufsort/).

## Usage instructions

The _gt_sain_plain_sortsuffixes_ function can be used to return a pointer to an integer array containing the offsets of the sorted suffixes. 

An example illustrating how to use it can be seen in example file [sk-sain-main.cpp](SAIS-SK/src/sk-sain-main.cpp). A [sample make file](SAIS-SK/src/Makefile) showing the necessary compilation and linkages is also provided.

``` bash
$ cd /tmp
$ git clone git@github.com:Ichindar/dstsc.git
$ cd dstsc/SAIS-SK/src/
$ make
$ ./sk-sain.x (some file)
```

Lempel-Ziv factorisor
---------------------

## Description 

Produces a Lempel-Ziv (LZ) factorisation of a byte string. A LZ factorisation is typically used for compression. It represents a string as a set of components, where future components are references to repeating previous strings.

For example, the LZ factorisation of the string 'ABCAB' is:

    ("A", 0), ("B", 0), ("C", 0), (0, 2)

Which means add characters a, b and c, and then go back to the first character and copy two characters, a and b.

Instructions can be self referencing. For instance, the LZ factorisation of 'ABCABCABC' is:

    ("A", 0), ("B", 0), ("C, 0), (0, 6)

## Implementation details

Employs a Suffix Array to produce a LZ factorisation of a byte stream. This implementation uses the LZ factorisation algorithm by Ohlebusch and Cog (2011), and is largerly based on code Simon Cog kindly provided. The algorithm has worse case complexity of O(n).

## Usage instructions

The _lz_factorise_file_ procedure can be used generate two vector arrays containing the lz factorisation of the contents of a file. An [example](lzOG/src/main.cpp) of how to use this function has also been provided.

The _lz_factorise_ procedure can be used to generate two vector arrays containg the lz factorisation of an unsigned char array. An [example](lzOG/src/main2.cpp) illustrating how to use this functon has been provided.

The _lz_refs_ procedure can be used to generate two vector arrays containg the lz factorisation of an unsigned char array. This is similar to the _lz_factorise_ procedure except it returns only the reference factors. For example, when parsing "ABCAB", the _lz_factorise_ procedure will return:

    (0, 2)

A [sample make file](lzOG/Makefile) has been provided showing the necessary compilation and linkage options required. 

``` bash
$ cd /tmp
$ git clone git@github.com:Ichindar/dstsc.git
$ cd dsts/lzOG/
$ make
$ bin/lz_factor src/hello
```

References
---------- 
  
  Enno Ohlebusch and Simon Cog (2011). Lempel-Ziv Factorization Revisited.
  Combinatorial Pattern Matching 2011 (CPM 2011). Lecture Notes for Computer Science.
  Springer Berlin Heidelberg. Volume 6661. pp 15-26.

  Ge Nong, Sem Zhang and Wai Hong Chan (2011). Two Efficient Algorithms for Linear Suffix Array Construction. _IEEE Transactions on Computers_. 60(10). pp 1471 - 1484.
