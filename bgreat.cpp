/*****************************************************************************
 *   Bgreat : De Bruijn Graph Read Mapping Tool
 *   Copyright (C) 2014  INRIA
 *   Authors: Antoine Limasset
 *   Contact: antoine.limasset@inria.fr, INRIA/IRISA/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
 *   Source: https://github.com/Malfoy/BGREAT
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
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



int main(int argc, char ** argv){
	// initRc();
	string reads, pairedReads, unitigs("unitig.fa"),pathFile("paths"), notAlignedFile("notAligned.fa");
	int errors(2), threads(1), ka(31), c, effort(2),dogMode(1);
	int anchorSize(ka);
	bool brute(false),fastq(false),correctionMode(false),orderKeep(false),vectorMode(false);
	while ((c = getopt (argc, argv, "u:x:k:g:m:t:e:f:a:i:bqcO")) != -1){
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
		case 'c':
			correctionMode=true;
			break;
		case 'O':
			orderKeep=true;
			break;
		}
	}
	if(not vectorMode){
		anchorSize=ka;
	}
	if(reads=="" and pairedReads==""){
		cout
		<<"-u read file (unpaired)"<<endl
		<<"-x read file (paired)"<<endl
		<<"-k k value (graph) (31)"<<endl
		<<"-a anchors length (k)"<<endl
		<<"-g unitig file (unitig.fa)"<<endl
		<<"-m number of missmatch allowed (2)"<<endl
		<<"-t number of threads (1)"<<endl
		<<"-e effort put in mapping (2)"<<endl
		<<"-f path file (paths)"<<endl
		<<"-q for fastq read file"<<endl
		<<"-c to output corrected reads"<<endl
		<<"-O to keep order of the reads"<<endl;
		return 0;
	}
	Aligner supervisor(unitigs,pathFile,notAlignedFile,ka,threads,errors,fastq,correctionMode,effort,dogMode,vectorMode,true,orderKeep,anchorSize);
	supervisor.indexUnitigs();
	if(reads!=""){
		supervisor.alignAll(not brute,reads,false);
	}
	if(pairedReads!=""){
		supervisor.alignAll(not brute,pairedReads,true);
	}
	return 0;
}
