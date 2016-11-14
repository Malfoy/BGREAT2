Warning: In developpement


BGREAT 2

bgeat now index kmers of the unitigs and do not need an other tool to align on large unitigs.
If too much memory is used try to reduce the proportion of kmer index with -i option -i 10 will index 1/10 kmers.
Bgreat index by default all kmers.
Bgreat can be used to correct reads from a DBG (correction mode) or to know where the reads appear in the graph.

usage:

-u read file (unpaired)
-x read file (paired)
-k k value (30)
-g unitig file (unitig.fa)
-m number of missmatch allowed (2)
-t number of thread (1)
-e effort put in mapping (2)
-f path file (paths)
-a not aligned file (notAligned.fa)
-q for fastq read file
-c to output corrected reads

