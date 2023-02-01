# ClaPIM
ClaPIM: Scalable Sequence CLAssification using Processing-In-Memory

This simulator aims to evaluate the improvment of adding a filter stage before the searching stage.
It checks the number of crossbar arrays that should have been checked and compare it to the number of crossbar arrays that will be checked after using the filter.
This reflects the improvment that will have on the power and the lifetime of the design.

To run the simulator you can run the current command:
> python filtering_improvment.py <file_name> <kmer_len> <kmers_inMAT> <eth> <num_queries>
  
filename:  *.fna file containing the genome
  
kmer_len:  the length of the kmers to extract from the genome
  
kmers_inMAT: number of kmers to fit in each crossbar array
  
eth:  the filter error threshold
  
num_queries: how many random queries to check. a larger number provides more precise average
