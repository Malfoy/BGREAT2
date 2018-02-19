# BGREAT2

## Improved version of BGREAT

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

[![Build Status](https://travis-ci.org/Malfoy/BWISE.svg?branch=master)](https://travis-ci.org/Malfoy/BGREAT2)


Bgreat now index anchors from the unitigs and do not need an other tool to align on large unitigs.
If too much memory is used try to reduce the proportion of kmer index with -i option as example -i 10 will index 1/10 kmers.
Bgreat index by default all kmers.
Bgreat can be used to correct reads from a DBG (correction mode) or to know where the reads appear in the graph.

## Usage:

-u read file
Unpaired Mapping

-x read file
Paired  mapping, reads should be interleaved

-k k value
Value of k used to construct the graph

-a anchors length
Size of the anchors used to start mapping
Can be used if k is way larger than 31

-g unitig file
Unitig file in fasta
This file can be obtained using bcalm (https://github.com/GATB/bcalm) on a reads file

-m number of missmatch allowed
Maximal hamming distance between a read and its corresponding graph sequence for a mapping to be valid
Default value is 5

-t number of thread

-f output file name

-q input reads are Fastq
Note that Bgreat ignore quality information

-c to output corrected reads
Bgreat will output the sequence of  corresponding path of the read in the graph
Intuitevely, the read is  "corrected" according to the graph sequence, mode used by Bcool corrector (https://github.com/Malfoy/BCOOL)

-O to keep read ordering

The advanced options are experimental and in current developpement and should not be used

## Path mode:

In the  default mode, the numbers outputed correspond to the paths of unitigs a read (or pair of reads) maps on.
`
>read1
3;4;-6; 
`

mean that the read1 mapped on unitig 3 then 4 then the reverse complement of the unitig 6.

To get the corresponding sequence the tool numberToSequences will do the conversion (warning: large files may be produced this way due to redundancy of large unitigs)

Usage:
./numbersToSequences  unitigs.fa paths 31 > superReads.fa



## In correction mode:
In this mode the corrected reads are direclty outputed.
If the -O option is used, the corrected reads will be in the right order.


## Example command lines:

Map an unpaired reads file on a low k DBG in a output file "output_paths"

./bgreat -u reads.fa  -g dbg27.fa -k 27 -f output_paths

Map an unpaired reads file on a low k DBG in a output file "output_paths" with a maximum of 2 missmatches

./bgreat -u reads.fa  -g dbg27.fa -k 27 -f output_paths -m 2

Map an unpaired reads file in FASTQ  on a low k DBG in a output file "output_paths"

./bgreat -u reads.fa -q -g dbg27.fa -k 27 -f output_paths


Map an unpaired reads file on a low k DBG in a output file "output_paths" using 8 cores

./bgreat -u reads.fa  -g dbg27.fa -k 27 -f output_paths -t 8


Map a   paired reads file (interleaved format) on a low k DBG in a output file "output_paths"

./bgreat -x paired_reads.fa  -g dbg27.fa -k 27 -f output_paths


Map a paired reads file (interleaved format) on a high k DBG in a output file "output_paths"  with a anchors size of 31 (good value for NGS reads)

./bgreat -x paired_reads.fa  -g dbg91.fa -k 91  -f output_paths -a 31



Correct an unpaired reads file on a low k DBG in a output file "output_paths"

./bgreat -u reads.fa  -g dbg27.fa -k 27 -f output_paths -c


Create superReads from a paired reads file on a low k DBG in a output file "output_paths"

./bgreat -x paired_reads.fa  -g dbg27.fa -k 27 -f output_paths -c














