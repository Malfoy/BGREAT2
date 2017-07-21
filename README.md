Improved version of BGREAT

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

[![Build Status](https://travis-ci.org/Malfoy/BWISE.svg?branch=master)](https://travis-ci.org/Malfoy/BGREAT2)




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
-q for fastq read file
-c to output corrected reads
-O to keep read ordering

In correction mode, only unpaired reads can be used, and corrected reads are outputed, if the -O option is used, the corrected reads will be in the right order.
In default mode, the numbers outputed correspond to the paths of unitigs the reads or pair of reads maps on.
>read1
3;4;-6; mean that the read1 mapped on unitig 3 then 4 then the reverse complement of the unitig 6.

To get the corresponding sequence the tool numberToSequences will do the conversion (warning: large file may be produced this way due to redundancy of large unitigs)

usage:
./numbersToSequences  unitigs.fa paths 31 > superReads.fa

