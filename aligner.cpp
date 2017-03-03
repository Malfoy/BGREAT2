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
 *  but WITHOUT ANY WARRANTY; without even the implied warranty ofF
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/



#include "aligner.h"
#include "utils.h"
#include "alignerExhaustive.cpp"
#include "alignerGreedy.cpp"
#include <thread>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>



using namespace std;



//Parse reads
void Aligner::getReads(vector<pair<string,string>>& reads, uint n){
	reads={};
	string read,header,inter;
	char c;
	if(fastq){
		for(uint i(0);i<n;++i){
			getline(readFile,header);
			getline(readFile,read);
			if(read.size()>2){
				bool fail(false);
				for(uint j(0);(j)<read.size();++j){
					if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G' and read[j]!='N'){
						fail=true;
						break;
					}
				}
				if(!fail){
					reads.push_back({header,read});
				}
			}
			getline(readFile,header);
			getline(readFile,header);
			if(readFile.eof()){return;}
		}
	}else{
		for(uint i(0);i<n;++i){
			getline(readFile,header);
			getline(readFile,read);
		point:
			c=readFile.peek();
			if(c=='>'){
				if(read.size()>2){
					bool fail(false);
					for(uint j(0);(j)<read.size();++j){
						if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G' and read[j]!='N'){
							fail=true;
							break;
						}
					}
					if(!fail){
						//~ if(read.size()>k){
							reads.push_back({header,read});
						//~ }
					}
				}
				read="";
			}else{
				if(!readFile.eof()){
					getline(readFile,inter);
					read+=inter;
					goto point;
				}else{
					if(read.size()>2){
						bool fail(false);
						for(uint j(0);(j)<read.size();++j){
							if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G' and read[j]!='N'){
								fail=true;
								break;
							}
						}
						if(!fail){
							//~ if(read.size()>k){
								reads.push_back({header,read});
							//~ }
						}
					}
					return;
				}
			}
		}
	}
}



kmer Aligner::getRepresentNum(const string& str){
	// string rc(reverseComplement(str));
	kmer a(str2num(str));
	kmer b(rcb(a,k-1));
	return ((b<=a) ? b : a);
}



string Aligner::num2str(kmer num){
	string str;
	int nuc;
	for(uint i(0);i<k-1;++i){
		nuc=num%4;
		switch (nuc){
			case 0:str.push_back('A');break;
			case 1:str.push_back('C');break;
			case 2:str.push_back('G');break;
			case 3:str.push_back('T');break;
		}
		num>>=2;
	}
	reverse( str.begin(), str.begin());
	return str;
}



//TODO replace vector by struct
vector<pair<string,uNumber>> Aligner::getEnd(kmer bin){
	vector<pair<string,uNumber>> result;
	kmer rc(rcb(bin,k-1));
	string unitig;
	unitigIndices indices;
	uNumber num;
	bool go(false);
	if(bin<rc){
		uint64_t hash=rightMPHF.lookup(bin);
		if(hash!=ULLONG_MAX){
			indices=rightIndices[hash];
			if(indices.overlap==bin){
				go=true;
			}
		}
	}else{
		uint64_t hash=leftMPHF.lookup(rc);
		if(hash!=ULLONG_MAX){
			indices=leftIndices[hash];
			if(indices.overlap==rc){
				go=true;
			}
		}
	}
	if(go){
		if(indices.indice1!=0){
			unitig=unitigs[indices.indice1];
			if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
				result.push_back({unitig,indices.indice1});
			}else{
				if(rcMode){
					result.push_back({unitigsRC[indices.indice1],-indices.indice1});
				}else{
					result.push_back({reverseComplements(unitig),-indices.indice1});
				}
			}
			if(indices.indice2!=0){
				unitig=unitigs[indices.indice2];
				if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
					result.push_back({unitig,indices.indice2});
				}else{
					if(rcMode){
						result.push_back({unitigsRC[indices.indice2],-indices.indice2});
					}else{
						result.push_back({reverseComplements(unitig),-indices.indice2});
					}
				}
				if(indices.indice3!=0){
					unitig=unitigs[indices.indice3];
					if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
						result.push_back({unitig,indices.indice3});
					}else{
						if(rcMode){
							result.push_back({unitigsRC[indices.indice3],-indices.indice3});
						}else{
							result.push_back({reverseComplements(unitig),-indices.indice3});
						}
					}
					if(indices.indice4!=0){
						unitig=unitigs[indices.indice4];
						if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
							result.push_back({unitig,indices.indice4});
						}else{
							if(rcMode){
								result.push_back({unitigsRC[indices.indice4],-indices.indice4});
							}else{
								result.push_back({reverseComplements(unitig),-indices.indice4});
							}
						}
					}
				}
			}
		}
	}
	return result;
}



vector<pair<string,uNumber>> Aligner::getEnd(string bin){
	vector<pair<string,uNumber>> result;
	string rc(reverseComplements(bin));
	string unitig;
	unitigIndicesstr indices;
	uNumber num;
	bool go(false);
	if(bin<rc){
		uint64_t hash=rightMPHFstr.lookup(bin);
		if(hash!=ULLONG_MAX){
			indices=rightIndicesstr[hash];
			if(indices.overlap==bin){
				go=true;
			}
		}
	}else{
		uint64_t hash=leftMPHFstr.lookup(rc);
		if(hash!=ULLONG_MAX){
			indices=leftIndicesstr[hash];
			if(indices.overlap==rc){
				go=true;
			}
		}
	}
	if(go){
		if(indices.indice1!=0){
			unitig=unitigs[indices.indice1];
			if((unitig.substr(unitig.size()-k+1,k-1))==bin){
				result.push_back({unitig,indices.indice1});
			}else{
				if(rcMode){
					result.push_back({unitigsRC[indices.indice1],-indices.indice1});
				}else{
					result.push_back({reverseComplements(unitig),-indices.indice1});
				}
			}
			if(indices.indice2!=0){
				unitig=unitigs[indices.indice2];
				if((unitig.substr(unitig.size()-k+1,k-1))==bin){
					result.push_back({unitig,indices.indice2});
				}else{
					if(rcMode){
						result.push_back({unitigsRC[indices.indice2],-indices.indice2});
					}else{
						result.push_back({reverseComplements(unitig),-indices.indice2});
					}
				}
				if(indices.indice3!=0){
					unitig=unitigs[indices.indice3];
					if((unitig.substr(unitig.size()-k+1,k-1))==bin){
						result.push_back({unitig,indices.indice3});
					}else{
						if(rcMode){
							result.push_back({unitigsRC[indices.indice3],-indices.indice3});
						}else{
							result.push_back({reverseComplements(unitig),-indices.indice3});
						}
					}
					if(indices.indice4!=0){
						unitig=unitigs[indices.indice4];
						if((unitig.substr(unitig.size()-k+1,k-1))==bin){
							result.push_back({unitig,indices.indice4});
						}else{
							if(rcMode){
								result.push_back({unitigsRC[indices.indice4],-indices.indice4});
							}else{
								result.push_back({reverseComplements(unitig),-indices.indice4});
							}
						}
					}
				}
			}
		}
	}
	return result;
}



vector<pair<string,uNumber>> Aligner::getBegin(kmer bin){
	vector<pair<string,uNumber>> result;
	kmer rc(rcb(bin,k-1));
	string unitig;
	unitigIndices indices;
	bool go(false);
	if(bin<rc){
		uint64_t hash=leftMPHF.lookup(bin);
		if(hash!=ULLONG_MAX){
			indices=leftIndices[hash];
			if(indices.overlap==bin){
				go=true;
			}
		}
	}else{
		uint64_t hash=rightMPHF.lookup(rc);
		if(hash!=ULLONG_MAX){
			indices=rightIndices[hash];
			if(indices.overlap==rc){
				go=true;
			}
		}
	}
	if(go){
		if(indices.indice1!=0){
			unitig=unitigs[indices.indice1];
			if(str2num(unitig.substr(0,k-1))==bin){
				result.push_back({unitig,indices.indice1});
			}else{
				if(rcMode){
					result.push_back({unitigsRC[indices.indice1],-indices.indice1});
				}else{
					result.push_back({reverseComplements(unitig),-indices.indice1});
				}
			}
			if(indices.indice2!=0){
				unitig=unitigs[indices.indice2];
				if(str2num(unitig.substr(0,k-1))==bin){
					result.push_back({unitig,indices.indice2});
				}else{
					if(rcMode){
						result.push_back({unitigsRC[indices.indice2],-indices.indice2});
					}else{
						result.push_back({reverseComplements(unitig),-indices.indice2});
					}
				}
				if(indices.indice3!=0){
					unitig=unitigs[indices.indice3];
					if(str2num(unitig.substr(0,k-1))==bin){
						result.push_back({unitig,indices.indice3});
					}else{
						if(rcMode){
							result.push_back({unitigsRC[indices.indice3],-indices.indice3});
						}else{
							result.push_back({reverseComplements(unitig),-indices.indice3});
						}
					}
					if(indices.indice4!=0){
						unitig=unitigs[indices.indice4];
						if(str2num(unitig.substr(0,k-1))==bin){
							result.push_back({unitig,indices.indice4});
						}else{
							if(rcMode){
								result.push_back({unitigsRC[indices.indice4],-indices.indice4});
							}else{
								result.push_back({reverseComplements(unitig),-indices.indice4});
							}
						}
					}
				}
			}
		}
	}
	return result;
}



vector<pair<string,uNumber>> Aligner::getBegin(string bin){
	vector<pair<string,uNumber>> result;
	string rc(reverseComplements(bin));
	string unitig;
	unitigIndicesstr indices;
	bool go(false);
	if(bin<rc){
		uint64_t hash=leftMPHFstr.lookup(bin);
		if(hash!=ULLONG_MAX){
			indices=leftIndicesstr[hash];
			if(indices.overlap==bin){
				go=true;
			}
		}
	}else{
		uint64_t hash=rightMPHFstr.lookup(rc);
		if(hash!=ULLONG_MAX){
			indices=rightIndicesstr[hash];
			if(indices.overlap==rc){
				go=true;
			}
		}
	}
	if(go){
		if(indices.indice1!=0){
			unitig=unitigs[indices.indice1];
			if((unitig.substr(0,k-1))==bin){
				result.push_back({unitig,indices.indice1});
			}else{
				if(rcMode){
					result.push_back({unitigsRC[indices.indice1],-indices.indice1});
				}else{
					result.push_back({reverseComplements(unitig),-indices.indice1});
				}
			}
			if(indices.indice2!=0){
				unitig=unitigs[indices.indice2];
				if((unitig.substr(0,k-1))==bin){
					result.push_back({unitig,indices.indice2});
				}else{
					if(rcMode){
						result.push_back({unitigsRC[indices.indice2],-indices.indice2});
					}else{
						result.push_back({reverseComplements(unitig),-indices.indice2});
					}
				}
				if(indices.indice3!=0){
					unitig=unitigs[indices.indice3];
					if((unitig.substr(0,k-1))==bin){
						result.push_back({unitig,indices.indice3});
					}else{
						if(rcMode){
							result.push_back({unitigsRC[indices.indice3],-indices.indice3});
						}else{
							result.push_back({reverseComplements(unitig),-indices.indice3});
						}
					}
					if(indices.indice4!=0){
						unitig=unitigs[indices.indice4];
						if((unitig.substr(0,k-1))==bin){
							result.push_back({unitig,indices.indice4});
						}else{
							if(rcMode){
								result.push_back({unitigsRC[indices.indice4],-indices.indice4});
							}else{
								result.push_back({reverseComplements(unitig),-indices.indice4});
							}
						}
					}
				}
			}
		}
	}
	return result;
}



string Aligner::recoverPath(vector<uNumber>& numbers,uint size){
	int offset(numbers[0]);
	string path(getUnitig(numbers[1]));
	for(uint i(2); i<numbers.size(); ++i){
		string unitig(getUnitig(numbers[i])),inter(compactionEndNoRC(path, unitig, k-1));
		if(inter.empty()){
			cout<<"bug compaction"<<endl;
			return {};
		}else{
			path=inter;
		}
	}
	path=path.substr(offset,size);
	return path;
}



string Aligner::recoverSuperReads(const vector<uNumber>& numbers){
	if(numbers.size()<1){
		return "";
	}
	string path(getUnitig(numbers[0])),unitig,inter;
	for(uint i(1); i<numbers.size(); ++i){
		unitig=(getUnitig(numbers[i]));inter=(compactionEndNoRC(path, unitig, k-1));
		if(inter.empty()){
			cout<<i<<endl;
			cout<<"bug compaction super reads"<<endl;
			return {};
		}else{
			path=inter;
		}
	}
	return path;
}



string Aligner::recoverSuperReadsNoStr(const vector<uNumber>& numbers){
	string path;
	if(numbers.size()<1){
		return "";
	}
	for(uint i(0); i<numbers.size(); ++i){
		path+=to_string(numbers[i])+";";
	}
	return path;
}



vector<uNumber> getcleanPaths(const vector<uNumber>& numbers, bool reverse,bool clean){
	vector<uNumber> res;
	if(clean){
		res=vector<uNumber>(&numbers[1],&numbers[numbers.size()]);
	}else{
		res=vector<uNumber>(&numbers[0],&numbers[numbers.size()]);
	}
	if(reverse){
		for(uint i(0); i<res.size();++i){
			res[i]=-res[i];
		}
		uNumber temp;
		for(uint i(0); i<res.size()/2;++i){
			temp=res[i];
			res[i]=res[res.size()-i-1];
			res[res.size()-i-1]=temp;
		}
	}
	return res;
}



bool Aligner::isNeighboor(const uint number1, const uint number2){
	string unitig1(getUnitig(number1)), unitig2(getUnitig(number2)),inter(compactionEndNoRC(unitig1, unitig2, k-1));
	if(inter.empty()){
		return false;
	}
	return true;
}



pair<string,string> Aligner::recoverSuperReadsPaired( const vector<uNumber>& vec, vector<uNumber>& vec2){
	//If one is empty return one read
	if(vec.size()<=0){
		if(vec2.size()<=0){
			return {"",""};
		}else{
			vector<uNumber> numbers2(getcleanPaths(vec2,false,true));
			return{recoverSuperReads(numbers2),""};
		}
	}else{
		if(vec2.size()<=0){
			vector<uNumber> numbers(getcleanPaths(vec,false,true));
			return{recoverSuperReads(numbers),""};
		}
	}

	//if they overlap
	vector<uNumber> numbers(getcleanPaths(vec,false,true));
	vector<uNumber> numbers2(getcleanPaths(vec2,true,true));
	for(uint i(0);i<numbers.size();++i){
		bool overlap(true);
		uint j(0);
		for(;j+i<numbers.size() and j<numbers2.size();++j){
			if(numbers[i+j]!=numbers2[j]){
				overlap=false;
				break;
			}
		}
		if(overlap){
			numbers.insert(numbers.end(),numbers2.begin()+j,numbers2.end());
			++superReads;
			return{recoverSuperReads(numbers),""};
		}
	}

	//it they do not overlap but can be compacted
	if(isNeighboor(numbers[numbers.size()-1],numbers2[0])){
		numbers.insert(numbers.end(),numbers2.begin(),numbers2.end());
		++superReads;
		return{recoverSuperReads(numbers),""};
	}
	return{recoverSuperReads(numbers),recoverSuperReads(numbers2)};
}



bool Aligner::compactVectors( vector<uNumber>& numbers, vector<uNumber>& numbers2){
	for(uint i(0);i<numbers.size();++i){
		bool overlap(true);
		uint j(0);
		for(;j+i<numbers.size() and j<numbers2.size();++j){
			if(numbers[i+j]!=numbers2[j]){
				overlap=false;
				break;
			}
		}
		if(overlap){
			numbers.insert(numbers.end(),numbers2.begin()+j,numbers2.end());
			numbers2={};
			return true;
		}
	}
	if(isNeighboor(numbers[numbers.size()-1],numbers2[0])){
		numbers.insert(numbers.end(),numbers2.begin(),numbers2.end());
		numbers2={};
			return true;
	}
}



pair<string,string> Aligner::recoverSuperReadsPairedNoStr( const vector<uNumber>& vec, vector<uNumber>& vec2){
	//If one is empty return one read
	if(vec.size()<=0){
		if(vec2.size()<=0){
			return {"",""};
		}else{
			vector<uNumber> numbers2(getcleanPaths(vec2,false,true));
			return{recoverSuperReadsNoStr(numbers2),""};
		}
	}else{
		if(vec2.size()<=0){
			vector<uNumber> numbers(getcleanPaths(vec,false,true));
			return{recoverSuperReadsNoStr(numbers),""};
		}
	}

	vector<uNumber> numbers(getcleanPaths(vec,false,true));
	vector<uNumber> numbers2(getcleanPaths(vec2,true,true));
	//if they overlap
	if(compactVectors(numbers,numbers2)){
		++superReads;
		return{recoverSuperReadsNoStr(numbers),""};
	}

	return{recoverSuperReadsNoStr(numbers),recoverSuperReadsNoStr(numbers2)};
}



string Aligner::getUnitig(int position){
	if(fullMemory){
		if(position>0){
			return unitigs[position];
		}
		if(rcMode){
			return unitigsRC[-position];
		}
		return reverseComplements(unitigs[-position]);
	}
	return "";

}



void Aligner::update(kmer& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdateOverlap;
}



void Aligner::updateK(kmer& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdateAnchor;
}



void Aligner::updateRC(kmer& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*k-4));
}



void Aligner::updateRCK(kmer& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*anchorSize-2));
}



vector<pair<kmer,uint>> Aligner::getListOverlap(const string& read){
	vector<pair<kmer,uint>> listOverlap;
	string overlap(read.substr(0,k-1));
	kmer num(str2num(overlap)),rcnum(rcb(num,k-1)), rep(min(num, rcnum));
	uint64_t hash;
	for(uint i(0);;++i){
		hash=leftMPHF.lookup(rep);
		if(true){//TODO test if is a true overlap
			listOverlap.push_back({num,i});
		}else{
			hash=rightMPHF.lookup(rep);
			if(true){
				listOverlap.push_back({num,i});
			}
		}
		if(i+k-1<read.size()){
			update(num,read[i+k-1]);
			updateRC(rcnum,read[i+k-1]);
			rep=(min(num, rcnum));
		}else{
			return listOverlap;
		}
	}
	return listOverlap;
}



vector<pair<kmer,uint>> Aligner::getNOverlap(const string& read, uint n){
	vector<pair<kmer,uint>> listOverlap;
	uint64_t hash;
	kmer num(str2num(read.substr(0,k-1))),rcnum(rcb(num,k-1)), rep(min(num, rcnum));
	for(uint i(0);;++i){
		bool done(false);
		hash=leftMPHF.lookup(rep);
		if(hash!=ULLONG_MAX){
			if(leftIndices[hash].overlap==rep){
				listOverlap.push_back({num,i});
				done=true;
			}
		}
		if(not done){
			hash=rightMPHF.lookup(rep);
			if(hash!=ULLONG_MAX){
				if(rightIndices[hash].overlap==rep){
					listOverlap.push_back({num,i});
				}
			}
		}
		if(listOverlap.size()>=n){
			return listOverlap;
		}
		if(i+k-1<read.size()){
			update(num,read[i+k-1]);
			updateRC(rcnum,read[i+k-1]);
			rep=(min(num, rcnum));
		}else{
			return listOverlap;
		}
	}
	return listOverlap;
}



vector<pair<pair<uint,uint>,uint>> Aligner::getNAnchors(const string& read, uint n){
	if(stringModeAnchor){
		return getNAnchorsstr(read,n);
	}
	return getNAnchorsnostr(read,n);


}



vector<pair<pair<uint,uint>,uint>> Aligner::getNAnchorsnostr(const string& read, uint n){
	unordered_set<uint> unitigsSelected;
	vector<pair<pair<uint,uint>,uint>> list;
	uint64_t hash;
	string unitig;
	kmer num(0),rcnum(0),rep(0);
	for(uint i(0);i+anchorSize<read.size();++i){
		bool returned(false);
		if(num==0 and rcnum==0){
			num=(str2num(read.substr(i,anchorSize)));
			rcnum=(rcb(num,anchorSize));
			rep=(min(num, rcnum));
		}else{
			updateK(num,read[i+anchorSize-1]);
			updateRCK(rcnum,read[i+anchorSize-1]);
			rep=(min(num, rcnum));
		}
		hash=anchorsMPHF.lookup(rep);
		if(hash!=ULLONG_MAX and anchorsChecking[hash]==rep){
			if(vectorMode){
				if(num==rep){
					for(uint j(0);j<anchorsPositionVector[hash].size();++j){
						if(unitigsSelected.count(anchorsPositionVector[hash][j].first)==0){
							unitigsSelected.insert(anchorsPositionVector[hash][j].first);
							list.push_back({anchorsPositionVector[hash][j],i});
						}
					}
				}else{
					for(uint j(0);j<anchorsPositionVector[hash].size();++j){
						if(unitigsSelected.count(anchorsPositionVector[hash][j].first)==0){
							unitigsSelected.insert(anchorsPositionVector[hash][j].first);
							list.push_back({{-anchorsPositionVector[hash][j].first,anchorsPositionVector[hash][j].second},i});
						}
					}
				}
			}else{
				if(unitigsSelected.count(anchorsPosition[hash].first)==0){
					unitigsSelected.insert(anchorsPosition[hash].first);
					if(num==rep){
						list.push_back({anchorsPosition[hash],i});
					}else{
						list.push_back({{-anchorsPosition[hash].first,anchorsPosition[hash].second},i});
					}
				}
			}
		}
		if(list.size()>=n){
			return list;
		}
	}
	return list;
}



vector<pair<pair<uint,uint>,uint>> Aligner::getNAnchorsstr(const string& read,uint n){
	unordered_set<uint> unitigsSelected;
	vector<pair<pair<uint,uint>,uint>> list;
	uint64_t hash;
	string unitig,num,rcnum,rep;
	uint positionUnitig;
	for(uint i(0);i+anchorSize<read.size();++i){
		bool returned(false);
		if(num.empty()){
			num=((read.substr(i,anchorSize)));
			rcnum=(reverseComplements(num));
			rep=(min(num, rcnum));
		}else{
			num=num.substr(1,anchorSize-1)+read[i+anchorSize];
			rcnum=revCompChar(read[i+anchorSize])+rcnum.substr(0,anchorSize-1);
			rep=(min(num,rcnum));
		}
		hash=anchorsMPHFstr.lookup(rep);
		if(hash!=ULLONG_MAX){
			if(vectorMode){
				if(num==rep){
					for(uint j(0);j<anchorsPositionVector[hash].size();++j){
						if(unitigsSelected.count(anchorsPositionVector[hash][j].first)==0){
							unitigsSelected.insert(anchorsPositionVector[hash][j].first);
							list.push_back({anchorsPositionVector[hash][j],i});
						}
					}
				}else{
					for(uint j(0);j<anchorsPositionVector[hash].size();++j){
						if(unitigsSelected.count(anchorsPositionVector[hash][j].first)==0){
							unitigsSelected.insert(anchorsPositionVector[hash][j].first);
							list.push_back({{-anchorsPositionVector[hash][j].first,anchorsPositionVector[hash][j].second},i});
						}
					}
				}
			}else{
				if(unitigsSelected.count(anchorsPosition[hash].first)==0){
					unitigsSelected.insert(anchorsPosition[hash].first);
					if(anchorsCheckingstr[hash]==rep){
						if(num==rep){
							list.push_back({anchorsPosition[hash],i});
						}else{
							list.push_back({{-anchorsPosition[hash].first,anchorsPosition[hash].second},i});
						}
					}
				}
			}
		}
		if(list.size()>=n){
			//~ cout<<"end"<<endl;
			return list;
		}
	}
	return list;
}


//TODO MULTITHREAD
void Aligner::indexUnitigsAux(){
	string line;
	unitigs.push_back("");
	unitigsRC.push_back("");
	uint leftsize,rightsize,anchorNumber;
	vector<kmer>* leftOver=new vector<kmer>;
	vector<kmer>* rightOver=new vector<kmer>;
	vector<kmer>* anchors=new vector<kmer>;
	while(!unitigFile.eof()){
		getline(unitigFile,line);
		getline(unitigFile,line);
		if(line.size()<k){
			break;
		}else{
			unitigs.push_back(line);
			unitigsRC.push_back(reverseComplements(line));
			kmer beg(str2num(line.substr(0,k-1))),rcBeg(rcb(beg,k-1));
			if(beg<=rcBeg){
				leftOver->push_back(beg);
			}else{
				rightOver->push_back(rcBeg);
			}
			kmer end(str2num(line.substr(line.size()-k+1,k-1))),rcEnd(rcb(end,k-1));
			if(end<=rcEnd){
				rightOver->push_back(end);
			}else{
				leftOver->push_back(rcEnd);
			}
			if(dogMode){
				kmer seq(str2num(line.substr(0,anchorSize))),rcSeq(rcb(seq,anchorSize)),canon(min(seq,rcSeq));
				anchors->push_back(canon);
				for(uint j(0);j+anchorSize<line.size();++j){
					updateK(seq,line[j+anchorSize]);
					updateRCK(rcSeq,line[j+anchorSize]);
					canon=(min(seq, rcSeq));
					anchors->push_back(canon);
				}
			}
		}
	}
	sort( leftOver->begin(), leftOver->end() );
	leftOver->erase( unique( leftOver->begin(), leftOver->end() ), leftOver->end() );
	sort( rightOver->begin(), rightOver->end() );
	rightOver->erase( unique( rightOver->begin(), rightOver->end() ), rightOver->end() );
	auto data_iterator = boomphf::range(static_cast<const kmer*>(&((*leftOver)[0])), static_cast<const kmer*>((&(*leftOver)[0])+leftOver->size()));
	leftMPHF= boomphf::mphf<kmer,hasher>(leftOver->size(),data_iterator,coreNumber,gammaFactor,false);
	leftsize=leftOver->size();
	delete leftOver;
	auto data_iterator2 = boomphf::range(static_cast<const kmer*>(&(*rightOver)[0]), static_cast<const kmer*>((&(*rightOver)[0])+rightOver->size()));
	rightMPHF= boomphf::mphf<kmer,hasher>(rightOver->size(),data_iterator2,coreNumber,gammaFactor,false);
	rightsize=rightOver->size();
	delete rightOver;
	if(dogMode){
		if(vectorMode){
			//TODO MULTITHREAD
			sort( anchors->begin(), anchors->end() );
			anchors->erase( unique( anchors->begin(), anchors->end() ), anchors->end() );
		}
		auto data_iterator3 = boomphf::range(static_cast<const kmer*>(&(*anchors)[0]), static_cast<const kmer*>((&(*anchors)[0])+anchors->size()));
		anchorsMPHF= boomphf::mphf<kmer,hasher>(anchors->size(),data_iterator3,coreNumber,gammaFactor,false);
	}
	anchorNumber=anchors->size();
	delete anchors;
	if(vectorMode){
		anchorsPositionVector.resize(anchorNumber,{});
	}else{
		anchorsPosition.resize(anchorNumber,{0,0});
	}
	anchorsChecking.resize(anchorNumber,0);
	leftIndices.resize(leftsize,{});
	rightIndices.resize(rightsize,{});
	fillIndices();
}


//TODO MULTITHREAD
void Aligner::indexUnitigsAuxStrfull(){
	string line,beg,end,seq,rcBeg,rcEnd,rcSeq,canon;
	unitigs.push_back("");
	unitigsRC.push_back("");
	uint leftsize,rightsize,anchorNumber;
	vector<string>* leftOver=new vector<string>;
	vector<string>* rightOver=new vector<string>;
	vector<string>* anchors=new vector<string>;
	while(!unitigFile.eof()){
		getline(unitigFile,line);
		getline(unitigFile,line);
		if(line.size()<k){
			break;
		}else{
			unitigs.push_back(line);
			unitigsRC.push_back(reverseComplements(line));
			beg=((line.substr(0,k-1)));
			rcBeg=(reverseComplements(beg));
			if(beg<=rcBeg){
				leftOver->push_back(beg);
			}else{
				rightOver->push_back(rcBeg);
			}
			end=(line.substr(line.size()-k+1,k-1));
			rcEnd=(reverseComplements(end));
			if(end<=rcEnd){
				rightOver->push_back(end);
			}else{
				leftOver->push_back(rcEnd);
			}
			if(dogMode){
				seq=((line.substr(0,anchorSize)));
				rcSeq=(reverseComplements(seq));
				canon=(min(seq,rcSeq));
				anchors->push_back(canon);
				for(uint j(0);j+anchorSize<line.size();++j){
					seq=seq.substr(1,anchorSize-1)+line[j+anchorSize];
					rcSeq=revCompChar(line[j+anchorSize])+rcSeq.substr(0,anchorSize-1);
					canon=(min(seq,rcSeq));
					anchors->push_back(canon);
				}
			}
		}
	}
	sort( leftOver->begin(), leftOver->end() );
	leftOver->erase( unique( leftOver->begin(), leftOver->end() ), leftOver->end() );
	sort( rightOver->begin(), rightOver->end() );
	rightOver->erase( unique( rightOver->begin(), rightOver->end() ), rightOver->end() );
	auto data_iterator = boomphf::range(static_cast<const string*>(&((*leftOver)[0])), static_cast<const string*>((&(*leftOver)[0])+leftOver->size()));
	leftMPHFstr= MPHFSTR(leftOver->size(),data_iterator,coreNumber,gammaFactor,false);
	leftsize=leftOver->size();
	delete leftOver;
	auto data_iterator2 = boomphf::range(static_cast<const string*>(&(*rightOver)[0]), static_cast<const string*>((&(*rightOver)[0])+rightOver->size()));
	rightMPHFstr= MPHFSTR(rightOver->size(),data_iterator2,coreNumber,gammaFactor,false);
	rightsize=rightOver->size();
	delete rightOver;
	if(dogMode){
		if(vectorMode){
			//TODO MULTITHREAD
			sort( anchors->begin(), anchors->end() );
			anchors->erase( unique( anchors->begin(), anchors->end() ), anchors->end() );
		}
		auto data_iterator3 = boomphf::range(static_cast<const string*>(&(*anchors)[0]), static_cast<const string*>((&(*anchors)[0])+anchors->size()));
		anchorsMPHFstr= MPHFSTR(anchors->size(),data_iterator3,coreNumber,gammaFactor,false);
	}
	anchorNumber=anchors->size();
	delete anchors;
	if(vectorMode){
		anchorsPositionVector.resize(anchorNumber,{});
	}else{
		anchorsPosition.resize(anchorNumber,{0,0});
	}
	leftIndicesstr.resize(leftsize,{});
	rightIndicesstr.resize(rightsize,{});
	anchorsCheckingstr.resize(anchorNumber,"");
	fillIndicesstr();
}



void Aligner::indexUnitigsAuxStrbutanchors(){
	string line,beg,end,seq,rcBeg,rcEnd,rcSeq,canon;
	unitigs.push_back("");
	unitigsRC.push_back("");
	uint leftsize,rightsize,anchorNumber;
	vector<string>* leftOver=new vector<string>;
	vector<string>* rightOver=new vector<string>;
	vector<kmer>* anchors=new vector<kmer>;
	while(!unitigFile.eof()){
		getline(unitigFile,line);
		getline(unitigFile,line);
		if(line.size()<k){
			break;
		}else{
			unitigs.push_back(line);
			unitigsRC.push_back(reverseComplements(line));
			beg=((line.substr(0,k-1)));
			rcBeg=(reverseComplements(beg));
			if(beg<=rcBeg){
				leftOver->push_back(beg);
			}else{
				rightOver->push_back(rcBeg);
			}
			end=(line.substr(line.size()-k+1,k-1));
			rcEnd=(reverseComplements(end));
			if(end<=rcEnd){
				rightOver->push_back(end);
			}else{
				leftOver->push_back(rcEnd);
			}
			if(dogMode){
				kmer seq(str2num(line.substr(0,anchorSize))),rcSeq(rcb(seq,anchorSize)),canon(min(seq,rcSeq));
				anchors->push_back(canon);
				for(uint j(0);j+anchorSize<line.size();++j){
					updateK(seq,line[j+anchorSize]);
					updateRCK(rcSeq,line[j+anchorSize]);
					canon=(min(seq, rcSeq));
					anchors->push_back(canon);
				}
			}
		}
	}
	sort( leftOver->begin(), leftOver->end() );
	leftOver->erase( unique( leftOver->begin(), leftOver->end() ), leftOver->end() );
	sort( rightOver->begin(), rightOver->end() );
	rightOver->erase( unique( rightOver->begin(), rightOver->end() ), rightOver->end() );
	auto data_iterator = boomphf::range(static_cast<const string*>(&((*leftOver)[0])), static_cast<const string*>((&(*leftOver)[0])+leftOver->size()));
	leftMPHFstr= MPHFSTR(leftOver->size(),data_iterator,coreNumber,gammaFactor,false);
	leftsize=leftOver->size();
	delete leftOver;
	auto data_iterator2 = boomphf::range(static_cast<const string*>(&(*rightOver)[0]), static_cast<const string*>((&(*rightOver)[0])+rightOver->size()));
	rightMPHFstr= MPHFSTR(rightOver->size(),data_iterator2,coreNumber,gammaFactor,false);
	rightsize=rightOver->size();
	delete rightOver;
	if(dogMode){
		if(vectorMode){
			//TODO MULTITHREAD
			sort( anchors->begin(), anchors->end() );
			anchors->erase( unique( anchors->begin(), anchors->end() ), anchors->end() );
		}
		auto data_iterator3 = boomphf::range(static_cast<const kmer*>(&(*anchors)[0]), static_cast<const kmer*>((&(*anchors)[0])+anchors->size()));
		anchorsMPHF= MPHF(anchors->size(),data_iterator3,coreNumber,gammaFactor,false);
	}
	anchorNumber=anchors->size();
	delete anchors;
	if(vectorMode){
		anchorsPositionVector.resize(anchorNumber,{});
	}else{
		anchorsPosition.resize(anchorNumber,{0,0});
	}
	leftIndicesstr.resize(leftsize,{});
	rightIndicesstr.resize(rightsize,{});
	anchorsChecking.resize(anchorNumber,0);
	fillIndicesstrbutanchors();
}



//TODO multihread
void Aligner::fillIndices(){
	unitigIndices indices;
	uint i;
	#pragma omp parallel for num_threads(coreNumber)
	for(i=1;i<unitigs.size();++i){
		string line(unitigs[i]);
		if(dogMode){
			kmer seq(str2num(line.substr(0,anchorSize))),rcSeq(rcb(seq,anchorSize)),canon(min(seq,rcSeq));
			uint64_t hash=anchorsMPHF.lookup(canon);
			if(canon==seq){
				if(vectorMode){
					mutexV[hash%1000].lock();
					anchorsPositionVector[hash].push_back({i,0});
					mutexV[hash%1000].unlock();
				}else{
					anchorsPosition[hash]={i,0};
				}
			}else{
				if(vectorMode){
					mutexV[hash%1000].lock();
					anchorsPositionVector[hash].push_back({-i,0});
					mutexV[hash%1000].unlock();
				}else{
					anchorsPosition[hash]={-i,0};
				}
			}
			anchorsChecking[hash]=canon;
			for(uint j(0);j+anchorSize<line.size();++j){
				updateK(seq,line[j+anchorSize]);
				updateRCK(rcSeq,line[j+anchorSize]);
				canon=(min(seq, rcSeq));
				uint64_t hash=anchorsMPHF.lookup(canon);
				if(canon==seq){
					if(vectorMode){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({i,j+1});
						mutexV[hash%1000].unlock();
					}else{
						anchorsPosition[hash]={i,j+1};
					}
				}else{
					if(vectorMode){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({-i,j+1});
						mutexV[hash%1000].unlock();
					}else{
						anchorsPosition[hash]={-i,j+1};
					}
				}
				anchorsChecking[hash]=canon;
			}
		}
	}
	string line;
	for(i=1;i<unitigs.size();++i){
		line=unitigs[i];
		kmer beg(str2num(line.substr(0,k-1))),rcBeg(rcb(beg,k-1));
		if(beg<=rcBeg){
			indices=leftIndices[leftMPHF.lookup(beg)];
			indices.overlap=beg;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			leftIndices[leftMPHF.lookup(beg)]=indices;
		}else{
			indices=rightIndices[rightMPHF.lookup(rcBeg)];
			indices.overlap=rcBeg;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			rightIndices[rightMPHF.lookup(rcBeg)]=indices;
		}
		kmer end(str2num(line.substr(line.size()-k+1,k-1))),rcEnd(rcb(end,k-1));
		if(end<=rcEnd){
			indices=rightIndices[rightMPHF.lookup(end)];
			indices.overlap=end;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			rightIndices[rightMPHF.lookup(end)]=indices;
		}else{
			indices=leftIndices[leftMPHF.lookup(rcEnd)];
			indices.overlap=rcEnd;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			leftIndices[leftMPHF.lookup(rcEnd)]=indices;
		}
	}
}



void Aligner::fillIndicesstr(){
	unitigIndicesstr indices;
	uint i;
	#pragma omp parallel for num_threads(coreNumber)
	for(i=(1);i<unitigs.size();++i){
		string line,seq,beg,rcBeg,rcEnd,end,rcSeq,canon;
		line=unitigs[i];
		if(dogMode){
			seq=((line.substr(0,anchorSize)));
			rcSeq=(reverseComplements(seq));
			canon=(min(seq,rcSeq));
			uint64_t hash=anchorsMPHFstr.lookup(canon);
			if(canon==seq){
				if(vectorMode){
					mutexV[hash%1000].lock();
					anchorsPositionVector[hash].push_back({i,0});
					mutexV[hash%1000].unlock();
				}else{
					anchorsPosition[hash]={i,0};
				}
			}else{
				if(vectorMode){
					mutexV[hash%1000].lock();
					anchorsPositionVector[hash].push_back({-i,0});
					mutexV[hash%1000].unlock();
				}else{
					anchorsPosition[hash]={-i,0};
				}
			}
			anchorsCheckingstr[hash]=canon;
			for(uint j(0);j+anchorSize<line.size();++j){
				seq=seq.substr(1,anchorSize-1)+line[j+anchorSize];
				rcSeq=revCompChar(line[j+anchorSize])+rcSeq.substr(0,anchorSize-1);
				canon=(min(seq,rcSeq));
				int64_t hash=anchorsMPHFstr.lookup(canon);
				if(canon==seq){
					if(vectorMode){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({i,j+1});
						mutexV[hash%1000].unlock();
					}else{
						anchorsPosition[hash]={i,j+1};
					}
				}else{
					if(vectorMode){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({-i,j+1});
						mutexV[hash%1000].unlock();
					}else{
						anchorsPosition[hash]={-i,j+1};
					}
				}
				anchorsCheckingstr[hash]=canon;
			}
		}
	}
	string line, seq,beg,rcBeg,rcEnd,end,rcSeq,canon;
	for(i=(1);i<unitigs.size();++i){
		line=unitigs[i];
		beg=((line.substr(0,k-1)));
		rcBeg=(reverseComplements(beg));
		if(beg<=rcBeg){
			indices=leftIndicesstr[leftMPHFstr.lookup(beg)];
			indices.overlap=beg;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			leftIndicesstr[leftMPHFstr.lookup(beg)]=indices;
		}else{
			indices=rightIndicesstr[rightMPHFstr.lookup(rcBeg)];
			indices.overlap=rcBeg;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			rightIndicesstr[rightMPHFstr.lookup(rcBeg)]=indices;
		}
		end=((line.substr(line.size()-k+1,k-1)));
		rcEnd=(reverseComplements(end));
		if(end<=rcEnd){
			indices=rightIndicesstr[rightMPHFstr.lookup(end)];
			indices.overlap=end;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			rightIndicesstr[rightMPHFstr.lookup(end)]=indices;
		}else{
			indices=leftIndicesstr[leftMPHFstr.lookup(rcEnd)];
			indices.overlap=rcEnd;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			leftIndicesstr[leftMPHFstr.lookup(rcEnd)]=indices;
		}
	}
}


void Aligner::fillIndicesstrbutanchors(){
	unitigIndicesstr indices;
	uint i;
	#pragma omp parallel for num_threads(coreNumber)
	for(i=(1);i<unitigs.size();++i){
		string line,seq,beg,rcBeg,rcEnd,end,rcSeq,canon;
		line=unitigs[i];
		if(dogMode){
			kmer seq(str2num(line.substr(0,anchorSize))),rcSeq(rcb(seq,anchorSize)),canon(min(seq,rcSeq));
			uint64_t hash=anchorsMPHF.lookup(canon);
			if(canon==seq){
				if(vectorMode){
					mutexV[hash%1000].lock();
					anchorsPositionVector[hash].push_back({i,0});
					mutexV[hash%1000].unlock();
				}else{
					anchorsPosition[hash]={i,0};
				}
			}else{
				if(vectorMode){
					mutexV[hash%1000].lock();
					anchorsPositionVector[hash].push_back({-i,0});
					mutexV[hash%1000].unlock();
				}else{
					anchorsPosition[hash]={-i,0};
				}
			}
			anchorsChecking[hash]=canon;
			for(uint j(0);j+anchorSize<line.size();++j){
				updateK(seq,line[j+anchorSize]);
				updateRCK(rcSeq,line[j+anchorSize]);
				canon=(min(seq, rcSeq));
				uint64_t hash=anchorsMPHF.lookup(canon);
				if(canon==seq){
					if(vectorMode){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({i,j+1});
						mutexV[hash%1000].unlock();
					}else{
						anchorsPosition[hash]={i,j+1};
					}
				}else{
					if(vectorMode){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({-i,j+1});
						mutexV[hash%1000].unlock();
					}else{
						anchorsPosition[hash]={-i,j+1};
					}
				}
				anchorsChecking[hash]=canon;
			}
		}
	}
	string line, seq,beg,rcBeg,rcEnd,end,rcSeq,canon;
	for(i=(1);i<unitigs.size();++i){
		line=unitigs[i];
		beg=((line.substr(0,k-1)));
		rcBeg=(reverseComplements(beg));
		if(beg<=rcBeg){
			indices=leftIndicesstr[leftMPHFstr.lookup(beg)];
			indices.overlap=beg;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			leftIndicesstr[leftMPHFstr.lookup(beg)]=indices;
		}else{
			indices=rightIndicesstr[rightMPHFstr.lookup(rcBeg)];
			indices.overlap=rcBeg;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			rightIndicesstr[rightMPHFstr.lookup(rcBeg)]=indices;
		}
		end=((line.substr(line.size()-k+1,k-1)));
		rcEnd=(reverseComplements(end));
		if(end<=rcEnd){
			indices=rightIndicesstr[rightMPHFstr.lookup(end)];
			indices.overlap=end;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			rightIndicesstr[rightMPHFstr.lookup(end)]=indices;
		}else{
			indices=leftIndicesstr[leftMPHFstr.lookup(rcEnd)];
			indices.overlap=rcEnd;
			if(indices.indice1==0){
				indices.indice1=i;
			}else if(indices.indice2==0){
				indices.indice2=i;
			}else if(indices.indice3==0){
				indices.indice3=i;
			}else{
				indices.indice4=i;
			}
			leftIndicesstr[leftMPHFstr.lookup(rcEnd)]=indices;
		}
	}
}


//TODO THIS WOLE PART COULD BE FACTORIZED
void Aligner::indexUnitigs(){
	auto startChrono=chrono::system_clock::now();
	uint nbThreads(1);
	vector<thread> threads;
	for (uint i(0); i<nbThreads; ++i){
		if(stringMode){
			if(stringModeAnchor){
				threads.push_back(thread(&Aligner::indexUnitigsAuxStrfull,this));
			}else{
				threads.push_back(thread(&Aligner::indexUnitigsAuxStrbutanchors,this));
			}
		}else{
			threads.push_back(thread(&Aligner::indexUnitigsAux,this));
		}
	}
	for(auto &t : threads){t.join();}
	auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
	cout<<"Indexing in seconds : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<endl;
}



void Aligner::alignAll(bool greedy, const string& reads, bool boolPaired){
	pairedMode=boolPaired;
	auto startChrono=chrono::system_clock::now();
	uint last(0);
	string file;
	for(uint i(0);i<reads.size();++i){
		if(reads[i]==','){
			file=reads.substr(last,i-last);
			readFile.close();
			readFile.open(file);
			cout<<file<<endl;
			unsigned char nbThreads(coreNumber);
			vector<thread> threads;
			for (size_t i(0); i<nbThreads; ++i){
				if(greedy){
					threads.push_back(thread(&Aligner::alignPartGreedy,this,i));
				}else{
					threads.push_back(thread(&Aligner::alignPartExhaustive,this));
				}
			}
			for(auto &t : threads){t.join();}
			last=i+1;
		}
	}
	file=reads.substr(last);
	readFile.close();
	readFile.open(file);
	cout<<file<<endl;
	unsigned char nbThreads(coreNumber);
	vector<thread> threads;
	for (size_t i(0); i<nbThreads; ++i){
		if(greedy){
			threads.push_back(thread(&Aligner::alignPartGreedy,this,i));
		}else{
			threads.push_back(thread(&Aligner::alignPartExhaustive,this));
		}
	}
	for(auto &t : threads){t.join();}

	cout<<"The End"<<endl;
	cout<<"Reads : "<<readNumber<<endl;
	cout<<"Not anchored : "<<noOverlapRead<<" Percent : "<<(100*float(noOverlapRead))/readNumber<<endl;
	cout<<"Anchored and aligned : "<<alignedRead<<" Percent : "<<(100*float(alignedRead))/(alignedRead+notAligned)<<endl;
	cout<<"Anchored but not aligned : "<<notAligned<<" Percent : "<<(100*float(notAligned))/(alignedRead+notAligned)<<endl;
	auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
	cout<<"Reads/seconds : "<<readNumber/(chrono::duration_cast<chrono::seconds>(waitedFor).count()+1)<<endl;
	cout<<"Mapping in seconds : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<endl;
	cout<<"Super reads : "<<superReads<<endl;
}



string Aligner::printPath(const vector<int32_t>& path){
	string res;
	for(uint i(0); i<path.size(); ++i){
		// *file<<path[i]<<'.';
		res+=to_string(path[i]);
		res+='.';
	}
	res+='\n';
	return res;
}
