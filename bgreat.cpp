/*****************************************************************************
 *   Bgreat : De Bruijn Graph Read Mapping Tool
 *   Authors: Antoine Limasset
 *   Contact: antoine.limasset@gmail.com
 *   Source: https://github.com/Malfoy/BGREAT
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty ofF
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/



#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "aligner.h"



using namespace std;



//~ TODO 3 MODES
//~ 1) unique optimal mappings
//~ 2) optimal mappings
//~ 3) all mappings
int main(int argc, char ** argv){
	// initRc();
	string reads, pairedReads, unitigs("unitig.fa"),pathFile("paths"), notAlignedFile("notAligned.fa");
	int errors(5), threads(1), ka(31), c, effort(1000),dogMode(1);
	int anchorSize(ka),ocurence_anchors(8);
	bool brute(false),fastq(false),correctionMode(false),orderKeep(false),vectorMode(false),preciseOutput(false),multi(false),printAlignment(false),allOptimalMapping(false),allMapping(false),compressOutput(false),anyOptimalMapping(false),compressionMode(false),Bulles(false);
	float ratioe(0.5);
	while ((c = getopt (argc, argv, "u:x:k:g:m:t:e:f:a:i:r:o:bqcOpMPABCFzZ")) != -1){
	switch(c){
		case 'u':
			reads=optarg;
			break;
		case 'x':
			pairedReads=optarg;
			break;
		case 'k':
			ka=stoi(optarg);
			break;
		case 'g':
			unitigs=(optarg);
			break;
		case 'm':
			errors=stoi(optarg);
			break;
		case 't':
			threads=stoi(optarg);
			break;
		case 'e':
			effort=stoi(optarg);
			break;
		case 'r':
			ratioe=stof(optarg);
			break;
		case 'f':
			pathFile=(optarg);
			break;
		case 'a':
			anchorSize=stoi(optarg);
			vectorMode=true;
			break;
		case 'b':
			brute=(true);
			break;
		case 'q':
			fastq=(true);
			break;
		case 'i':
			dogMode=stoi(optarg);
			break;
		case 'o':
			ocurence_anchors=stoi(optarg);
			break;
		case 'c':
			correctionMode=true;
			break;
		case 'O':
			orderKeep=true;
			break;
		case 'p':
			preciseOutput=true;
			break;
		case 'M':
			multi=true;
			break;
		case 'P':
			printAlignment=true;
			break;
		case 'A':
			allMapping=true;
			break;
		case 'B':
			allOptimalMapping=true;
			break;
		case 'F':
			anyOptimalMapping=true;
			break;
		case 'C':
			compressionMode=true;
			break;
		case 'Z':
			Bulles=true;
			break;
		case 'z':
			compressOutput=true;
			break;
		}
	}
	if(not vectorMode){
		anchorSize=ka;
	}
	if(reads=="" and pairedReads==""){
		cout
		<<"Standard options"<<endl
		<<"-u read file (unpaired)"<<endl
		<<"-x read file (paired)"<<endl
		<<"-k k value (graph) (31)"<<endl
		<<"-a anchors length (k)"<<endl
		<<"-g unitig file (unitig.fa)"<<endl
		<<"-m number of missmatch allowed (5)"<<endl
		<<"-t number of threads (1)"<<endl
		<<"-f output file (paths)"<<endl
		<<"-q for fastq read file"<<endl
		<<"-O to keep order of the reads"<<endl<<endl

		<<"Advanced options"<<endl
		<<"-c to output corrected reads"<<endl
		<<"-C to output compressed reads"<<endl
		<<"-p to more precise output"<<endl
		<<"-P to print the alignments"<<endl
		<<"-A to output all possible mapping"<<endl
		<<"-B to output all possible optimal mapping"<<endl
		<<"-o maximal occurence of an anchor (8)"<<endl
		<<"-e effort put in mapping (1000)"<<endl
		<<"-F to output any optimal mapping"<<endl;
		//~ <<"-C to output any optimal mapping"<<endl;
		return 0;
	}
	Aligner supervisor(unitigs,pathFile,notAlignedFile,ka,threads,errors,fastq,correctionMode,effort,dogMode,vectorMode,true,orderKeep,anchorSize,preciseOutput,multi,ratioe,allOptimalMapping,allMapping,printAlignment,compressOutput,anyOptimalMapping,compressionMode,ocurence_anchors);
	supervisor.indexUnitigs();

	if(Bulles){
		supervisor.graphFile.open("popped_dbg.fa");
		supervisor.Crush_bubbles();
		return 0;
	}
	if(reads!=""){
		supervisor.alignAll(not brute,reads,false);
	}
	if(pairedReads!=""){
		supervisor.alignAll(not brute,pairedReads,true);
	}
	return 0;
}
