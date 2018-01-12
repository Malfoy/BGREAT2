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



#ifndef __BGREAT__Aligner__
#define __BGREAT__Aligner__



#include <stdio.h>
#include <fstream>
#include <vector>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include "utils.h"
#include "BooPHF.h"
#include "zstr.hpp"



using namespace std;



class Custom_string_Hasher
{
public:
	// the class should have operator () with this signature :
	uint64_t operator ()   (std::string key, uint64_t seed=0) const
	{
		uint64_t hash  =  hash_fn(key);
		hash ^= seed;
		return hash;
	}
     std::hash<std::string> hash_fn;
};


//then tell BBhash to use this custom hash : (also appears below, line 104)
typedef boomphf::mphf<  std::string, Custom_string_Hasher  > MPHFSTR;




typedef boomphf::SingleHashFunctor<kmer>  hasher;
typedef boomphf::mphf<  kmer, hasher  > MPHF;



namespace std { template <> struct hash<__uint128_t> {
	typedef __uint128_t argument_type;
	typedef uint64_t result_type; uint64_t operator()(__uint128_t key) const { return transform_to_size_t(key); } };
}


struct unitigIndices{
	kmer overlap;
	uint32_t indice1;
	uint32_t indice2;
	uint32_t indice3;
	uint32_t indice4;
};



struct unitigIndicesstr{
	string overlap;
	uint32_t indice1;
	uint32_t indice2;
	uint32_t indice3;
	uint32_t indice4;
};



struct unitigIndicesVector{
	kmer overlap;
	vector<uNumber> vec;
};






class Aligner{
public:
	ifstream unitigFile;
	istream* readFile;
	ofstream pathFile, noOverlapFile, notMappedFile;
	FILE * pathFilef;
	FILE * readFileF;
	vector<FILE *> pathFileComp;
	ostream* pathCompressed;
	//~ FILE * notMappedFilef;
	//~ FILE * readFileF;
	MPHF leftMPHF,rightMPHF,anchorsMPHF;
	MPHFSTR leftMPHFstr,rightMPHFstr,anchorsMPHFstr;
	vector<unitigIndices> leftIndices,rightIndices;
	vector<unitigIndicesstr> leftIndicesstr,rightIndicesstr;
	vector<pair<int32_t,uint32_t>> anchorsPosition;
	vector<vector<pair<int32_t,uint32_t>>> anchorsPositionVector;
	vector<uint8_t> anchorsChecking;
	//~ vector<string> anchorsCheckingstr;
	atomic<uint> alignedRead, readNumber, noOverlapRead, notAligned, unitigNumber, overlaps, iter, superReads,notCompatedSR, overlappingPath,overlappingStr,singleMiddle,neighbor,includedPath,failed_pair;
	vector<string> unitigs, unitigsRC;
	kmer offsetUpdateOverlap;
	kmer offsetUpdateAnchor;
	uint coreNumber, gammaFactor, errorsMax, tryNumber, fracKmer,k,threadToPrint,threadToRead,anchorSize,maxPositionAnchors, iterLoop, nbBuckets,sizeBuckets;
	double ratioError;
	mutex unitigMutex, unitigMutex2, readMutex, indexMutex, pathMutex, noOverlapMutex, notMappedMutex, progressMutex;
	array<mutex,1000> mutexV;
	array<mutex,1024> pathMutexComp;

	string unitigFileName, pathToWrite;
	bool correctionMode, vectorMode, rcMode, fastq, dogMode,fullMemory,pairedMode,stringMode,keepOrder, preciseOutput,stringModeAnchor,noMultiMapping,uniqueOptimalMappingMode,optimalMappingMode, printAlignment,headerNeeded,compression, anyOptimalMapping,compressionMode;

	Aligner(const string& Unitigs, const string& paths, const string& notMapped, uint kValue, uint cores,uint errorsAllowed, bool bfastq, bool bcorrectionMode, uint effort, uint dogModeInt, bool vectorModeBool, bool rcModeBool,bool orderKeep,uint anchorsSize,bool preciseB,bool multi,float ratioe,bool ballOptimalMapping,bool ballMapping,bool bprintAlignment,bool compressOutput,bool any,bool compressionM,uint occurence_anchors){
		nbBuckets=2;
		compressionMode=compressionM;
		iterLoop=0;
		if(any){
			noMultiMapping=false;
		}else{
			noMultiMapping=true;
		}
		maxPositionAnchors=occurence_anchors;
		compression=compressOutput;
		headerNeeded=false;
		uniqueOptimalMappingMode=true;
		optimalMappingMode=false;

		if(compression){
			pathCompressed= (new zstr::ofstream((paths+".gz").c_str()));
		}else{
			pathFilef=fopen(paths.c_str(),"wb");
		}
		if(ballOptimalMapping){
			uniqueOptimalMappingMode=false;
			optimalMappingMode=true;
		}
		if(ballMapping){
			uniqueOptimalMappingMode=false;
			optimalMappingMode=false;
		}

		if(compressionMode){
			for(uint i(0);i<nbBuckets;++i){
				pathFileComp.push_back(fopen((paths+to_string(i)).c_str(),"wb"));
			}
		}
		printAlignment=bprintAlignment;
		preciseOutput=preciseB;
		anchorSize=anchorsSize;
		keepOrder=orderKeep;
		unitigFileName=Unitigs;
		rcMode=rcModeBool;
		dogMode=fullMemory=true;
		unitigFile.open(unitigFileName);


		ratioError=ratioe;
		k=kValue;
		if(k>63){
			stringMode=true;
		}else{
			stringMode=false;
		}
		if(anchorsSize>=k){
			anchorSize=k;
			vectorMode=false;
		}else{
			vectorMode=true;
		}
		if(anchorSize>63){
			stringModeAnchor=true;
		}else{
			stringModeAnchor=false;
		}
		coreNumber=cores;
		errorsMax=errorsAllowed;
		tryNumber=effort;
		gammaFactor=10;
		fracKmer=dogModeInt;
		correctionMode=bcorrectionMode;
		fastq=bfastq;
		notCompatedSR=threadToPrint=superReads=alignedRead=readNumber=noOverlapRead=notAligned=unitigNumber=overlaps=0;
		offsetUpdateOverlap=1;
		offsetUpdateOverlap<<=(2*(k-1));
		offsetUpdateAnchor=1;
		offsetUpdateAnchor<<=(2*(anchorSize));
		iter=1;
	}

	~Aligner(){
		if(compression){
			delete(pathCompressed);
		}else{
			fclose(pathFilef);
		}
		unitigFile.close();
	}

	void indexUnitigs();
	void alignAll(bool, const string&,bool);
	void alignPartGreedy(uint indice);
	void alignPartExhaustive();
	void indexUnitigsAux();
	vector<pair<kmer,uint>> getListOverlap(const string& read);
	uint checkEndExhaustive(const string& read, const pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	bool cover(const string& read, vector<pair<kmer,uint>>& listOverlap, uint start, uint end,vector<uNumber>& path);
	bool cover2(const string& read, const vector<pair<kmer,uint>>& listOverlap, const uint start, uint end,vector<uNumber>& path);
	uint mapOnRightAll(const string &read, vector<uNumber>& path, pair<kmer, uint> overlap, const vector<pair<kmer,uint>>& listOverlap, bool& end, uint start);
	bool mapOnRightEndAll(const string &read, vector<uNumber>& path, pair<kmer, uint> overlap);
	bool mapOnRightEndMax(const string &read, vector<uNumber>& path, pair<kmer, uint> overlap);
	bool mapOnLeftEndAll(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap);
	bool mapOnLeftEndFirst(const string &read, vector<uNumber>& path,const pair<kmer, uint>& overlap);
	bool mapOnLeftEndMax(const string &read, vector<uNumber>& path, const  pair<kmer, uint>& overlap);
	uint mapOnLeftEndExhaustive(const string &read, vector<uNumber>& path, const pair<kmer, uint>& , uint errors);
	uint mapOnRightEndExhaustive(const string &read, vector<uNumber>& path, const pair<kmer, uint>& , uint errors);
	string getUnitig(int position);
	string getUnitigFile(uint position);
	void getReads(vector<pair<string,string>>& reads, uint n);
	vector<uNumber> alignReadGreedy(const string& read, bool& overlapFound, uint errors,bool& rc);
	uint checkBeginExhaustive(const string& read, const pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	uint checkPair(const pair<kmer, uint>& overlap1, const pair<kmer, uint>& overlap2, const string& read,uNumber& path,uint errorsAllowed);
	uint coverGreedy(const string& read, const vector<pair<kmer,uint>>& listOverlap, const uint start, uint end, vector<uNumber>& path, uint  errors, bool&);
	pair<uint,int> mapOnRight2(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,uint start, uint errors);
	uint mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	uint checkBeginGreedy(const string& read,const pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	uint mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	vector<uNumber> alignReadExhaustive(const string& read, bool& overlapFound, uint errors);
	uint checkEndGreedy(const string& read,const pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	string readUnitig(uint position);
	vector<pair<string,uNumber>> getBegin(kmer);
	vector<pair<string,uNumber>> getEnd(kmer);
	string num2str(kmer num);
	kmer getRepresentNum(const string& str);
	string recoverPath(vector<uNumber>& numbers,uint size);
	vector<uNumber> getBeginNumber(kmer bin);
	vector<uNumber> getEndNumber(kmer bin);
	void updateRC(kmer&	min, char nuc);
	void update(kmer&	min, char nuc);
	void updateRCK(kmer&	min, char nuc);
	void updateK(kmer&	min, char nuc);
	pair<uint,uint> mapOnRight(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,uint start, uint errors);
	uint mapOnRightEndExhaustivePartial(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	uint mapOnLeftEndExhaustivePartial(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	uint coverGreedyPath(const string& read, const vector<pair<kmer,uint>>& listOverlap, const uint start, uint end, vector<uNumber>& path, uint  errors, bool& ended);
	uint checkPairPaths(const pair<kmer, uint>& overlap1, const pair<kmer, uint>& overlap2, const string& read, uNumber& number, uint errorsAllowed, const vector<uNumber> path);
	vector<uNumber> alignReadGreedyPath(const string& read, bool& overlapFound, uint errors, bool rc);
	uint checkBeginExhaustivePath(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	uint checkEndExhaustivePath(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	pair<uint,uint> mapOnRightPath(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,uint start, uint errors);
	uint mapOnLeftEndExhaustivePaths(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	uint mapOnRightEndExhaustivePath(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	vector<uNumber> alignReadExhaustivePath(const string& read, bool& overlapFound, uint errors);
	vector<overlapStruct> getListOverlapStuct(const string& read);
	vector<uNumber> alignReadGreedyCache(const string& read, bool& overlapFound, uint errors, bool rc);
	uint coverGreedyCache(const string& read, const vector<overlapStruct>& listOverlap, const uint start, uint end, vector<uNumber>& path, uint  errors, bool& ended);
	vector<overlapStruct> getListOverlapCache(const string& read);
	uint checkEndGreedyCache(const string& read, overlapStruct& overlap, vector<uNumber>& path, uint errors);
	uint checkBeginGreedyCache(const string& read, overlapStruct& overlap, vector<uNumber>& path, uint errors);
	uint checkPairCache(const overlapStruct& overlap1, const overlapStruct& overlap2, const string& read, uNumber& number, uint errorsAllowed);
	pair<uint,uint> mapOnRightCache(const string &read, vector<uNumber>& path, const overlapStruct& overlap, const  vector<overlapStruct>& listOverlap, bool& ended,uint start, uint errors);
	vector<pair<kmer,uint>> getNOverlap(const string& read, uint n);
	vector<uNumber> alignReadExhaustiveR(const string& read, bool& overlapFound, uint errors);
	vector<uNumber> alignReadGreedyCover(const string& read, bool& overlapFound, uint errors, bool rc);
	vector<string> getEndStr(kmer bin);
	vector<string> getBeginStr(kmer bin);
	void knowNeighbour();
	vector<pair<string,uNumber>> getBeginOpti(kmer bin, uNumber last);
	vector<pair<string,uNumber>> getEndOpti(kmer bin, uNumber last);
	string printPath(const vector<int32_t>& path);
	vector<uNumber> alignReadGreedyAnchors(const string& read, bool& overlapFound, uint errorMax, bool& rc, bool& nooverlap);
	vector<pair<pair<uint,uint>,uint>> getNAnchors(const string& read,uint n);
	string recoverSuperReads(const vector<uNumber>& numbers);
	pair<string,string> recoverSuperReadsPaired( const vector<uNumber>& vec, vector<uNumber>& vec2);
	string recoverSuperReadsNoStr(const vector<uNumber>& numbers, uint offset);
	//~ string recoverSuperReadsNoStr(const vector<uNumber>& numbers);
	pair<string,string> recoverSuperReadsPairedNoStr( const vector<uNumber>& vec, vector<uNumber>& vec2);
	bool isNeighboor(const uint number1, const uint number2);
	void fillIndices();
	void fillIndicesVector();
	vector<pair<string,uNumber>> getBeginV(kmer bin);
	vector<pair<string,uNumber>> getEndV(kmer bin);
	void indexUnitigsAuxStr();
	void fillIndicesstr();
	vector<pair<pair<uint,uint>,uint>> getNAnchorsstr(const string& read,uint n);
	vector<pair<string,uNumber>> getBegin(string bin);
	vector<pair<string,uNumber>> getEnd(string bin);
	vector<uNumber> alignReadGreedyAnchorsstr(const string& read, bool& overlapFound, uint errorMax, bool& rc, bool& noOverlap);
	uint checkBeginGreedy(const string& read,const pair<string, uint>& overlap, vector<uNumber>& path, uint errors);
	uint checkEndGreedy(const string& read,const pair<string, uint>& overlap, vector<uNumber>& path, uint errors);
	uint mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<string, uint>& overlap , uint errors);
	uint mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<string, uint>& overlap , uint errors);
	bool compactVectors( vector<uNumber>& numbers1, vector<uNumber>& numbers2);
	void indexUnitigsAuxStrbutanchors();
	void indexUnitigsAuxStrfull();
	void fillIndicesstrbutanchors();
	vector<pair<pair<uint,uint>,uint>> getNAnchorsnostr(const string& read, uint n);
	vector<uNumber> alignReadGreedyAnchors(const string& read, uint errorMax,const pair<pair<uint,uint>,uint>& anchor, uint& error);;
	vector<uNumber> alignReadGreedyAnchorsstr(const string& read, uint errorMax, const pair<pair<uint,uint>,uint>& anchor,uint& error);
	void alignReadOpti(const string& read, vector<int>& path,bool perfect);
	vector<int> inclued(vector<int>& v1, vector<int>& v2);
	void getReads2(vector<pair<string,string>>& reads, uint n);
	uint missmatchNumber(const string& seq1, const string& seq2, unsigned int n);
	void alignReadFrom(const string& read, vector<int>& path, int unumber);
	string recoverSuperReadsCor(const vector<uNumber>& numbers,uint readSize);
	void alignReadAllOpti(const string& read, vector<vector<int>>& pathVector);
	void alignReadAll(const string& read, vector<vector<int>>& pathVector);
	vector<uNumber> cleanSR(const vector<uNumber>& numbers, uint readSize);
	vector<uNumber> getcleanPaths(const vector<uNumber>& numbers, bool reverse,bool clean);
	void addIndicesEnd(unitigIndices& indices, kmer bin,vector<pair<string,uNumber>>& result);
	void addIndicesBegin(unitigIndices& indices, kmer bin,vector<pair<string,uNumber>>& result);
	void addIndicesEndStr(unitigIndicesstr& indices, string bin,vector<pair<string,uNumber>>& result);
	void addIndicesBeginStr(unitigIndicesstr& indices, string bin,vector<pair<string,uNumber>>& result);
	string path2nuc(const vector<int32_t>& path);
	vector<uNumber> path_clean(const vector<uNumber>& numbers, uint readSize);
	uint find_path_to(uNumber numbers, uNumber numbers2,vector<uNumber>& res, uint depth,int start);
	void Crush_bubbles();


};



#endif /* defined(__BGREAT__Aligner__) */
