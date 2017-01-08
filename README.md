A DNA K-mer is a substring (or a sequence) of length K

Given an uncompressed FASTQ file, this program produces a list of the most frequent K-mers, sorted by frequency for any K or list size.

"FASTQ files are used for storing both a biological sequence and its corresponding quality scores"
(https://en.wikipedia.org/wiki/FASTQ_format)

The files themselves can be extremely large.

Usage:

		./kmercounter --f [filename] --k [k-mer length] --l [output list size] --i [implementation]

a file name following option --f must be specified

default k-mer length is 30

default output list size is 25

default implementation is bf for Bloom Filter

other valid options are:

		sm: a two step set map based implementation
	
		mo: a 1 step map only implementation
  
general guidelines for choosing an impelentation:

		for any K larger than 15, choose bloom filter unless the input file size is very small (<=50mb)
		
		for any K smaller than 5, map only implementation is generally better even if the file size is very large
		
		in between above K values, two step set map implementation pulls ahead
  


