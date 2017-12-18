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



#include "aligner.h"
#include "utils.h"
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
using namespace chrono;



uint8_t getHash(const string& str){
	uint8_t result(0);
	for(uint i(0);i<4;++i){
		result+=nuc2int(str[i]);
	}
	return result;
}



uint8_t getHash(kmer n){
	return n%256;
}

//Parse reads
//TODO faster
void Aligner::getReads(vector<pair<string,string>>& reads, uint n){
	reads={};
	string read,header,inter;
	char c;
	if(fastq){
		for(uint i(0);i<n;++i){
			getline(*readFile,header);
			getline(*readFile,read);
			if(read.size()>0){
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
			getline(*readFile,header);
			getline(*readFile,header);
			if(readFile->eof()){return;}
		}
	}else{
		for(uint i(0);i<n;++i){
			getline(*readFile,header);
			getline(*readFile,read);
		point:
			c=readFile->peek();
			if(c=='>'){
				if(read.size()>0){
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
				if(!readFile->eof()){
					getline(*readFile,inter);
					read+=inter;
					goto point;
				}else{
					if(read.size()>0){
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


//TODO FA
void Aligner::getReads2(vector<pair<string,string>>& reads, uint n){
	reads={};
	string read,header,inter;
	int buffSize(1000000);
	char buff[buffSize];
	int realPred(-1),pred(-1);
	while (not feof(readFileF)){
		pred=-1;
		int res=fread(buff, 1, buffSize, readFileF);
		for(uint i(0);i<res;++i){
			if(buff[i]=='\n'){
				if(header==""){
					header=string(buff+pred+1,i-pred-1);
					pred=i;
				}else{
					read=string(buff+pred+1,i-pred-1);
					reads.push_back({header,read});
					header=read="";
					realPred=pred=i;
					if(reads.size()>=n){
						fseek(readFileF,realPred-res,SEEK_CUR);
						return;
					}
				}
			}else{
			}

		}
		if(-res+pred+1<0){
			fseek(readFileF,-res+pred+1,SEEK_CUR);
		}else{
		}
	}
}



uint Aligner::missmatchNumber(const string& seq1, const string& seq2, uint n){
	uint miss(0);
	for(uint i(0); i<seq2.size(); ++i){
		if(seq2[i]!=seq1[i]){
			++miss;
		}
	}
	if(((double)miss)/(double)seq2.size()>ratioError){
		return 1000;
	}
	return miss;
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
	//~ cout<<"ge"<<endl;
	vector<pair<string,uNumber>> result;
	kmer rc(rcb(bin,k-1));
	string unitig;
	unitigIndices indices;
	uNumber num;
	bool go(false);
	if(bin<=rc){
		uint64_t hash=rightMPHF.lookup(bin);
		if(hash!=ULLONG_MAX){
			indices=rightIndices[hash];
			if(indices.overlap==bin){
				addIndicesEnd(indices,bin,result);
			}
		}
	}
	if(bin>=rc){
		uint64_t hash=leftMPHF.lookup(rc);
		if(hash!=ULLONG_MAX){
			indices=leftIndices[hash];
			if(indices.overlap==rc){
				addIndicesEnd(indices,bin,result);
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
	if(bin<=rc){
		uint64_t hash=rightMPHFstr.lookup(bin);
		if(hash!=ULLONG_MAX){
			indices=rightIndicesstr[hash];
			if(indices.overlap==bin){
				addIndicesEndStr(indices,bin,result);
			}
		}
	}else{
		uint64_t hash=leftMPHFstr.lookup(rc);
		if(hash!=ULLONG_MAX){
			indices=leftIndicesstr[hash];
			if(indices.overlap==rc){
				addIndicesEndStr(indices,bin,result);
			}
		}
	}
	return result;
}



void Aligner::addIndicesBegin(unitigIndices& indices, kmer bin,vector<pair<string,uNumber>>& result){
	kmer rc(rcb(bin,k-1));
	if(indices.indice1!=0){
		string unitig=unitigs[indices.indice1];
		if(str2num(unitig.substr(0,k-1))==bin){
			result.push_back({unitig,indices.indice1});
		}
		if(str2num(unitig.substr(unitig.size()-k+1))==rc){
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
			}
			if(str2num(unitig.substr(unitig.size()-k+1))==rc){
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
				}
				if(str2num(unitig.substr(unitig.size()-k+1))==rc){
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
					}
					if(str2num(unitig.substr(unitig.size()-k+1))==rc){
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


void Aligner::addIndicesBeginStr(unitigIndicesstr& indices, string bin,vector<pair<string,uNumber>>& result){
	string rc(reverseComplements(bin));
	if(indices.indice1!=0){
		string unitig=unitigs[indices.indice1];
		if((unitig.substr(0,k-1))==bin){
			result.push_back({unitig,indices.indice1});
		}
		if((unitig.substr(unitig.size()-k+1))==rc){
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
			}
			if((unitig.substr(unitig.size()-k+1))==rc){
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
				}
				if((unitig.substr(unitig.size()-k+1))==rc){
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
					}
					if((unitig.substr(unitig.size()-k+1))==rc){
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


void Aligner::addIndicesEnd(unitigIndices& indices, kmer bin,vector<pair<string,uNumber>>& result){
	kmer rc(rcb(bin,k-1));
	if(indices.indice1!=0){
		string unitig=unitigs[indices.indice1];
		if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
			result.push_back({unitig,indices.indice1});
		}
		if(str2num(unitig.substr(0,k-1))==rc){
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
			}
			if(str2num(unitig.substr(0,k-1))==rc){
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
				}
				if(str2num(unitig.substr(0,k-1))==rc){
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
					}
					if(str2num(unitig.substr(0,k-1))==rc){
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



void Aligner::addIndicesEndStr(unitigIndicesstr& indices, string bin,vector<pair<string,uNumber>>& result){
	string rc(reverseComplements(bin));
	if(indices.indice1!=0){
		string unitig=unitigs[indices.indice1];
		if((unitig.substr(unitig.size()-k+1,k-1))==bin){
			result.push_back({unitig,indices.indice1});
		}
		if((unitig.substr(0,k-1))==rc){
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
			}
			if((unitig.substr(0,k-1))==rc){
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
				}
				if((unitig.substr(0,k-1))==rc){
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
					}
					if((unitig.substr(0,k-1))==rc){
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




vector<pair<string,uNumber>> Aligner::getBegin(kmer bin){
	vector<pair<string,uNumber>> result;
	kmer rc(rcb(bin,k-1));
	string unitig;
	unitigIndices indices;
	bool go(false);
	if(bin<=rc){
		uint64_t hash=leftMPHF.lookup(bin);
		if(hash!=ULLONG_MAX){
			indices=leftIndices[hash];
			if(indices.overlap==bin){
				addIndicesBegin(indices,bin,result);
			}
		}
	}
	if(bin>=rc){
		uint64_t hash=rightMPHF.lookup(rc);
		if(hash!=ULLONG_MAX){
			indices=rightIndices[hash];
			if(indices.overlap==rc){
				addIndicesBegin(indices,bin,result);
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
	if(bin<=rc){
		uint64_t hash=leftMPHFstr.lookup(bin);
		if(hash!=ULLONG_MAX){
			indices=leftIndicesstr[hash];
			if(indices.overlap==bin){
				addIndicesBeginStr(indices,bin,result);
			}
		}
	}
	if(bin>=rc){
		uint64_t hash=rightMPHFstr.lookup(rc);
		if(hash!=ULLONG_MAX){
			indices=rightIndicesstr[hash];
			if(indices.overlap==rc){
				addIndicesBeginStr(indices,bin,result);
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
			cin.get();
			return {};
		}else{
			path=inter;
		}
	}
	return path;
}



string Aligner::recoverSuperReadsCor(const vector<uNumber>& numbers, uint readSize){
	if(numbers.size()<2){
		return "";
	}
	string path(getUnitig(numbers[1])),unitig,inter;
	for(uint i(2); i<numbers.size(); ++i){
		unitig=(getUnitig(numbers[i]));inter=(compactionEndNoRC(path, unitig, k-1));
		if(inter.empty()){
			cout<<i<<endl;
			cout<<"bug compaction super reads"<<endl;
			return {};
		}else{
			path=inter;
		}
	}
	//~ cout<<"length: "<<path.size()<<endl;
	return path.substr(numbers[0], readSize);
}


vector<uNumber> Aligner::path_clean(const vector<uNumber>& numbers, uint readSize){
	vector<uNumber> res;
	if(numbers.size()<2 ){
		return res;
	}
	//~ cout<<"list"<<endl;
	//~ for(uint i(0);i<numbers.size();++i){
		//~ cout<<numbers[i]<<" ";
	//~ }
	//~ cout<<endl;
	uint pos_start(numbers[0]);
	uint indice_start(1);
	//~ cout<<"A"<<flush;
	while(pos_start+k-1 >= unitigs[abs(numbers[indice_start])].size() and indice_start+2<=numbers.size()){
		pos_start=k-1-(unitigs[abs(numbers[indice_start])].size()-pos_start);
		indice_start++;
	}
	//~ cout<<"B"<<flush;
	res.push_back(pos_start);
	res.push_back(numbers[indice_start]);
	string path(getUnitig(numbers[indice_start])),unitig,inter;
	//~ cout<<"D"<<endl;
	if(path.size()-pos_start>=readSize){
		return res;
	}
	//~ cout<<"C"<<flush;
	for(uint i(indice_start+1); i<numbers.size(); ++i){
		unitig=(getUnitig(numbers[i]));
		inter=(compactionEndNoRC(path, unitig, k-1));
		res.push_back(numbers[i]);
		if(inter.empty()){
			cout<<i<<endl;
			cout<<"bug compaction super reads"<<endl;
			cin.get();
			return {};
		}else{
			if(inter.size()-pos_start>=readSize){
				return res;
			}
			path=inter;
		}
	}
	//~ cout<<4<<flush;
	return res;
}



//TODO CHECK END ALSO
vector<uNumber> Aligner::cleanSR(const vector<uNumber>& numbers, uint readSize){
	//~ if(numbers.size()<3){return numbers;}
	vector<uNumber> result;
	string unitig;
	uint position(numbers[0]);
	uint lastPosition;
	uint i(1);
	for(;i<numbers.size();++i){
		unitig=(getUnitig(numbers[i]));
		if(position+k-1<unitig.size() or position+readSize<=unitig.size()){
			break;
		}else{
			lastPosition=position;
			position-=(unitig.size()-k+1);
		}
	}
	int lentgthAlignement(unitig.size()-position);
	result.push_back(position);
	result.push_back(numbers[i]);
	++i;
	for(;i<numbers.size();++i){
		if(lentgthAlignement<readSize){
			unitig=(getUnitig(numbers[i]));
			lentgthAlignement+(unitig.size()-k+1);
			result.push_back(numbers[i]);
		}else{
			break;
		}

	}

	string readCore1(recoverSuperReadsCor(numbers,readSize));
	string readCore2(recoverSuperReadsCor(result,readSize));
	//~ if(readCore1!=readCore2){
		//~ cout<<readCore1<<" and "<<readCore2<<endl;
	//~ }
	//~ if(result.size()<numbers.size()){
			//~ cout<<"win";
	//~ }

	return result;
}



string Aligner::recoverSuperReadsNoStr(const vector<uNumber>& numbers, uint offset=0){
	string path;
	if(numbers.size()<1+offset){
		return "";
	}
	for(uint i(offset); i<numbers.size(); ++i){
		path+=to_string(numbers[i])+";";
	}
	return path;
}



vector<uNumber> Aligner::getcleanPaths(const vector<uNumber>& numbers, bool reverse,bool clean){
	vector<uNumber> res;
	if(numbers.empty()){
		return res;
	}
	if(clean){
		if(numbers[0]+k-1> unitigs[abs(numbers[1])].size()){
			res=vector<uNumber>(&numbers[2],&numbers[numbers.size()]);
		}else{
			res=vector<uNumber>(&numbers[1],&numbers[numbers.size()]);
		}
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
	//~ cout<<"euuuuuh"<<endl;
	string unitig1(getUnitig(number1)), unitig2(getUnitig(number2));
	//~ cout<<"gocompaction"<<endl;
	string inter(compactionEndNoRC(unitig1, unitig2, k-1));
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



string overlapping(const string& str1, const string& str2, uint overlapMin){
	if(overlapMin>str1.size() or overlapMin>str2.size()){return "";}
	string suffix(str1.substr(str1.size()-overlapMin));
	int pos = str2.find(suffix, 0);
	if(pos !=-1){
		//~ int temp = str2.find(suffix,pos+1);
		//~ if(temp !=-1){
			//~ cout<<"fail1"<<endl;cin.get();
			//~ return "";
		//~ }%TODO PROBLEM OF SUFFIX REPETITION
	}else{
		return "";
	}
	if(overlapMin+pos<=str1.size()){

		if(str2.substr(0,pos)==str1.substr(str1.size()-overlapMin-pos,pos)){
			return str1+str2.substr((uint)pos+overlapMin);
		}
	}
	return "";
}



bool equalV(const vector<uNumber>& numbers,const vector<uNumber>& numbers2,int begin1, int begin2, int length){
	if(begin1<0 or begin2<0){
		return false;
	}
	for(uint i(0);i<length;++i){
		if(numbers[begin1+i]!=numbers[begin2+i]){
			return false;
		}
	}
	return true;
}



uint Aligner::find_path_to(uNumber numbers, uNumber numbers2, vector<uNumber>& res, uint depth, int max_size){
	//RETURN 0 if 0 PATH FOUND 1 IF ONE PAH FOUND 2 if MULTIPLE PATH FOUND
	string unitig(getUnitig(numbers));
	vector<pair<string,uNumber>> next;
	vector<uNumber> inter;
	vector<uNumber> recursion;
	vector<uNumber> res_sauv=res;
	bool found(false);
	uint valid;
	int next_unitig;
	if(stringMode){
		next=(getBegin((unitig.substr(unitig.size()-k+1))));
	}else{
		next=(getBegin(str2num(unitig.substr(unitig.size()-k+1))));
	}

	for(uint i(0);i<next.size();++i){
		if(next[i].second==numbers2){
			if(found){
				return 2;
			}else{
				found=true;
				inter=res;
			}
			//~ return 1;
			//~ cout<<endl;
			//~ for(uint iii(0);iii<res.size();++iii){
				//~ cout<<res[iii]<<" ";
			//~ }
			//~ cout<<endl;
		}else{
			if(depth<5){
				recursion.push_back(next[i].second);
			}
			//~ res=res_sauv;
			//~ res.push_back(next[i].second);
			//~ uint valid=find_path_to(next[i].second,numbers2,res,depth+1,start);
			//~ if (valid==2){
				//~ cout<<"ARGH"<<endl;
				//~ return 2;
			//~ }
			//~ if(valid==0){
				//~ res=res_sauv;
				//~ continue;
			//~ }
			//~ if(valid==1){

				//~ if(found){
					//~ cout<<"ARGH"<<endl;
					//~ return 2;
				//~ }else{
					//~ found=true;
					//~ inter=res;
				//~ }
			//~ }
		}
	}

	for(uint i(0);i<recursion.size();++i){
		res=res_sauv;
		res.push_back(recursion[i]);
		uint valid=0;
		if(unitigs[abs(recursion[i])].size()-k+1 < max_size){
			valid=find_path_to(recursion[i],numbers2,res,depth+1,max_size-(unitigs[abs(recursion[i])].size()-k+1));
		}
		if(valid==2){
			return 2;
		}
		if(valid==1){
			if(found){
				return 2;
			}else{
				found=true;
				inter=res;
			}
		}
	}

	if(found){
		res=inter;
		return 1;
	}
	return 0;
}




bool Aligner::compactVectors(vector<uNumber>& numbers, vector<uNumber>& numbers2){
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
			uint oldSize(numbers.size());
			uint oldSize2(numbers2.size());
			numbers.insert(numbers.end(),numbers2.begin()+j,numbers2.end());
			numbers2={};
			if(numbers.size()> oldSize and numbers.size()> oldSize2){
				overlappingPath++;
			}else{
				includedPath++;
			}
			return true;
		}
	}
	//~ string unitig(recoverSuperReads(numbers));
	//~ string unitig2((recoverSuperReads(numbers2)));
	//~ string merge(overlapping(unitig,unitig2,50));
	//~ if(merge!=""){
		//~ vector<uNumber> numbers3;
		//~ alignReadFrom(merge,numbers3,numbers[0]);
		//~ if(not numbers3.empty()){
			//~ numbers3=getcleanPaths(numbers3,false,true);
			//~ numbers=numbers3;
			//~ numbers2={};
			//~ overlappingStr++;

			//~ return true;
		//~ }
	//~ }

	vector<uNumber> inter;
	//~ cout<<"c parti"<<endl;
	uint res=find_path_to(numbers[numbers.size()-1],numbers2[0],inter,0,2000);

	//~ uint res=0;
	if(res==1){
		//~ cout<<"success"<<endl;

		//~ if(inter.size()>0){
			//~ cout<<numbers[numbers.size()-1]<<" "<<numbers2[0]<<endl;
			//~ cout<<"number ";
			//~ for(uint iii(0);iii<numbers.size();++iii){
				//~ cout<<numbers[iii]<<" ";
			//~ }
			//~ cout<<endl;
			//~ cout<<"number2 ";
			//~ for(uint iii(0);iii<numbers2.size();++iii){
				//~ cout<<numbers2[iii]<<" ";
			//~ }
			//~ cout<<endl;
			//~ cout<<"inter ";
			//~ for(uint iii(0);iii<inter.size();++iii){
				//~ cout<<inter[iii]<<" ";
			//~ }
			//~ cout<<endl;
			numbers.insert(numbers.end(),inter.begin(),inter.end());
			numbers.insert(numbers.end(),numbers2.begin(),numbers2.end());
			numbers2={};
			//~ for(uint iii(0);iii<numbers.size();++iii){
				//~ cout<<numbers[iii]<<" ";
			//~ }
			//~ cout<<endl<<endl;
			//~ recoverSuperReads(numbers);
			//~ if (numbers[numbers.size()-1]==18 or numbers[numbers.size()-1]==-18 or numbers2[0]==18 or numbers2[0]==-18){
				//~ cin.get();
			//~ }
		//~ }
		singleMiddle++;
		//~ cout<<3<<flush;
		return true;
	}
	//~ cout<<"c fail"<<endl;
	//~ cout<<"number ";
	//~ for(uint iii(0);iii<numbers.size();++iii){
		//~ cout<<numbers[iii]<<" ";
	//~ }
	//~ cout<<endl;
	//~ cout<<"number2 ";
	//~ for(uint iii(0);iii<numbers2.size();++iii){
		//~ cout<<numbers2[iii]<<" ";
	//~ }
	//~ cin.get();
	failed_pair++;
	return false;
}



pair<string,string> Aligner::recoverSuperReadsPairedNoStr( const vector<uNumber>& vec, vector<uNumber>& vec2){
	//If one is empty return one read
	vector<uNumber> numbers(getcleanPaths(vec,false,true));
	vector<uNumber> numbers2(getcleanPaths(vec2,true,true));
	//~ cout<<"clean"<<endl;
	if(numbers.size()<=0){
		if(numbers2.size()<=0){
			return {"",""};
		}else{
			return{recoverSuperReadsNoStr(numbers2),""};
		}
	}else{
		if(numbers2.size()<=0){
			return{recoverSuperReadsNoStr(numbers),""};
		}
	}
	//~ cout<<vec.size()<<endl;

	//if they overlap
	auto numberSAUV=numbers;
	if(compactVectors(numbers,numbers2)){
		++superReads;
		if("-18;-34;39;35;-18;38;"==recoverSuperReadsNoStr(numbers)){cout<<"found"<<endl;
			for(uint i(0);i<numberSAUV.size();++i){
				cout<<numberSAUV[i]<<";";
			}
			cout<<endl;
			for(uint i(0);i<numbers2.size();++i){
				cout<<numbers2[i]<<";";
			}
			cout<<"found"<<endl;
			exit(0);
		}
		return{recoverSuperReadsNoStr(numbers),""};
	}
	//~ cout<<"go"<<endl;
	//~ for(uint i(0);i<numbers.size();++i){
		//~ cout<<numbers[i]<<" ";
	//~ }
	//~ cout<<endl;
		//~ for(uint i(0);i<numbers2.size();++i){
		//~ cout<<numbers2[i]<<" ";
	//~ }
//~ cout<<endl;
	++notCompatedSR;
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
	unordered_map<uint,uint> unitigsSelected;
	//~ unordered_map<uint,uint> unitigsSelected2;
	vector<pair<pair<uint,uint>,uint>> list;
	uint64_t hash,anchorAdded(0);
	string unitig;
	kmer num(0),rcnum(0),rep(0);
	for(uint i(0); i+anchorSize<read.size(); ++i){
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
		if(hash!=ULLONG_MAX and anchorsChecking[hash]==getHash(rep)){
			bool addedAnAnchor(true);
			if(vectorMode){
				if(num==rep){
					for(uint j(0);j<anchorsPositionVector[hash].size();++j){
						int32_t unitigNum(anchorsPositionVector[hash][j].first);
						int32_t positionUnitig(anchorsPositionVector[hash][j].second);
						uint unitigNumPos(unitigNum>0?unitigNum:-unitigNum);
						if(positionUnitig<=unitigs[unitigNumPos].size()-k){
							if(unitigsSelected.count(unitigNumPos)==0){
								list.push_back({anchorsPositionVector[hash][j],i});
								//~ addedAnAnchor=true;
							}else{
								//~ if(unitigsSelected[unitigNumPos]!=positionUnitig+1 and unitigsSelected[unitigNumPos]!=positionUnitig-1){
									//~ list.push_back({anchorsPositionVector[hash][j],i});
								//~ }
							}
							unitigsSelected[unitigNumPos]=positionUnitig;
						}
					}
				}else{
					for(uint j(0);j<anchorsPositionVector[hash].size();++j){
						int32_t unitigNum(anchorsPositionVector[hash][j].first);
						int32_t positionUnitig(anchorsPositionVector[hash][j].second);
						uint unitigNumPos(unitigNum>0?unitigNum:-unitigNum);
						if(positionUnitig+anchorSize>=k){
							if(unitigsSelected.count(unitigNumPos)==0){
								list.push_back({{-anchorsPositionVector[hash][j].first,anchorsPositionVector[hash][j].second},i});
								//~ addedAnAnchor=true;
							}else{
								//~ if(unitigsSelected[unitigNumPos]!=positionUnitig+1 and unitigsSelected[unitigNumPos]!=positionUnitig-1){
									//~ list.push_back({{-anchorsPositionVector[hash][j].first,anchorsPositionVector[hash][j].second},i});
								//~ }
							}
							unitigsSelected[unitigNumPos]=positionUnitig;
						}
					}
				}
			}else{
				int32_t unitigNum(anchorsPosition[hash].first);
				uint unitigNumPos(unitigNum>0?unitigNum:-unitigNum);
				int32_t positionUnitig(anchorsPosition[hash].second);
				if(unitigsSelected.count(unitigNumPos)==0){
					if(num==rep){
							list.push_back({anchorsPosition[hash],i});
						}else{
							list.push_back({{-anchorsPosition[hash].first,anchorsPosition[hash].second},i});
						}
				}else{
					//~ if(unitigsSelected[unitigNumPos]!=positionUnitig+1 and unitigsSelected[unitigNumPos]!=positionUnitig-1){
						//~ if(num==rep){
							//~ list.push_back({anchorsPosition[hash],i});
						//~ }else{
							//~ list.push_back({{-anchorsPosition[hash].first,anchorsPosition[hash].second},i});
						//~ }
					//~ }
				}
				unitigsSelected[unitigNumPos]=positionUnitig;
			}
			if(addedAnAnchor){
				anchorAdded++;
			}
		}
		if(anchorAdded>=n){
			return list;
		}
	}
	return list;
}



vector<pair<pair<uint,uint>,uint>> Aligner::getNAnchorsstr(const string& read,uint n){
	unordered_map<uint,uint> unitigsSelected;
	vector<pair<pair<uint,uint>,uint>> list;
	uint64_t hash;
	string unitig,num,rcnum,rep;
	uint positionUnitig;
	uint anchorAdded(0);
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
		if(hash!=ULLONG_MAX and anchorsChecking[hash]==getHash(rep)){
			anchorAdded++;
			if(vectorMode){
				if(num==rep){
					for(uint j(0);j<anchorsPositionVector[hash].size();++j){
						int32_t unitigNum(anchorsPositionVector[hash][j].first);
						uint unitigNumPos(unitigNum>0?unitigNum:-unitigNum);
						if(unitigsSelected.count(unitigNumPos)==0){
							list.push_back({anchorsPositionVector[hash][j],i});
						}else{
							//~ if(unitigsSelected[unitigNumPos]!=positionUnitig+1 and unitigsSelected[unitigNumPos]!=positionUnitig-1){
								//~ list.push_back({anchorsPositionVector[hash][j],i});
							//~ }
						}
						unitigsSelected[unitigNumPos]=positionUnitig;
					}
				}else{
					for(uint j(0);j<anchorsPositionVector[hash].size();++j){
						int32_t unitigNum(anchorsPositionVector[hash][j].first);
						uint unitigNumPos(unitigNum>0?unitigNum:-unitigNum);
						if(unitigsSelected.count(unitigNumPos)==0){
							list.push_back({{-anchorsPositionVector[hash][j].first,anchorsPositionVector[hash][j].second},i});
						}else{
							//~ if(unitigsSelected[unitigNumPos]!=positionUnitig+1 and unitigsSelected[unitigNumPos]!=positionUnitig-1){
								//~ list.push_back({{-anchorsPositionVector[hash][j].first,anchorsPositionVector[hash][j].second},i});
							//~ }
						}
						unitigsSelected[unitigNumPos]=positionUnitig;
					}
				}
			}else{
				int32_t unitigNum(anchorsPosition[hash].first);
				uint unitigNumPos(unitigNum>0?unitigNum:-unitigNum);
				if(unitigsSelected.count(unitigNumPos)==0){
					unitigsSelected[unitigNumPos]=positionUnitig;
					if(num==rep){
						list.push_back({anchorsPosition[hash],i});
					}else{
						list.push_back({{-anchorsPosition[hash].first,anchorsPosition[hash].second},i});
					}
				}
			}
		}
		if(anchorAdded>=n){
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
	rightOver->push_back(0);
	leftOver->push_back(0);
	cout<<"Reading Unitigs: "<<flush;auto start1=system_clock::now();
	while(!unitigFile.eof()){
		getline(unitigFile,line);
		getline(unitigFile,line);
		if(line.size()<k){
			continue;
		}else{
			unitigs.push_back(line);
			unitigsRC.push_back(reverseComplements(line));
			kmer beg(str2num(line.substr(0,k-1))),rcBeg(rcb(beg,k-1));
			if(beg<=rcBeg){
				leftOver->push_back(beg);
			}
			if(beg>=rcBeg){
				rightOver->push_back(rcBeg);
			}
			kmer end(str2num(line.substr(line.size()-k+1,k-1))),rcEnd(rcb(end,k-1));
			if(end<=rcEnd){
				rightOver->push_back(end);
			}
			if(end>=rcEnd){
				leftOver->push_back(rcEnd);
			}
			if(dogMode){
				kmer seq(str2num(line.substr(0,anchorSize))),rcSeq(rcb(seq,anchorSize)),canon(min(seq,rcSeq));
				anchors->push_back(canon);
				for(uint j(0);j+anchorSize<line.size();++j){
					updateK(seq,line[j+anchorSize]);
					updateRCK(rcSeq,line[j+anchorSize]);
					if((j+1)%fracKmer==0){
						canon=(min(seq, rcSeq));
						anchors->push_back(canon);
					}
				}
			}
		}
	}
	auto end1=system_clock::now();auto waitedFor1=end1-start1;cout<<"Duration "<<duration_cast<seconds>(waitedFor1).count()<<" seconds"<<endl;

	cout<<"Sorting anchors: "<<flush;auto start2=system_clock::now();
	sort( leftOver->begin(), leftOver->end() );
	leftOver->erase( unique( leftOver->begin(), leftOver->end() ), leftOver->end() );
	sort( rightOver->begin(), rightOver->end() );
	rightOver->erase( unique( rightOver->begin(), rightOver->end() ), rightOver->end() );
	auto end2=system_clock::now();auto waitedFor2=end2-start2;cout<<"Duration "<<duration_cast<seconds>(waitedFor2).count()<<" seconds"<<endl;

	cout<<"Creating MPHF: "<<flush;auto start3=system_clock::now();
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
	auto end3=system_clock::now();auto waitedFor3=end3-start3;cout<<"Duration "<<duration_cast<seconds>(waitedFor3).count()<<" seconds"<<endl;

	cout<<"Filling index: "<<flush;auto start4=system_clock::now();
	fillIndices();
	auto end4=system_clock::now();auto waitedFor4=end4-start4;cout<<"Duration "<<duration_cast<seconds>(waitedFor4).count()<<" seconds"<<endl;
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
			continue;
		}else{
			unitigs.push_back(line);
			unitigsRC.push_back(reverseComplements(line));
			beg=((line.substr(0,k-1)));
			rcBeg=(reverseComplements(beg));
			if(beg<=rcBeg){
				leftOver->push_back(beg);
			}
			if(beg>=rcBeg){
				rightOver->push_back(rcBeg);
			}
			end=(line.substr(line.size()-k+1,k-1));
			rcEnd=(reverseComplements(end));
			if(end<=rcEnd){
				rightOver->push_back(end);
			}
			if(end>=rcEnd){
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
					if((j+1)%fracKmer==0){
						canon=(min(seq,rcSeq));
						anchors->push_back(canon);
					}
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
	anchorsChecking.resize(anchorNumber,0);
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
	leftOver->push_back("");
	rightOver->push_back("");
	cout<<"Reading Unitigs: "<<flush;auto start1=system_clock::now();
	while(!unitigFile.eof()){
		getline(unitigFile,line);
		getline(unitigFile,line);
		if(line.size()>=k){
			unitigs.push_back(line);
			unitigsRC.push_back(reverseComplements(line));
			beg=((line.substr(0,k-1)));
			rcBeg=(reverseComplements(beg));
			if(beg<=rcBeg){
				leftOver->push_back(beg);
			}
			if(beg>=rcBeg){
				rightOver->push_back(rcBeg);
			}
			end=(line.substr(line.size()-k+1,k-1));
			rcEnd=(reverseComplements(end));
			if(end<=rcEnd){
				rightOver->push_back(end);
			}
			if(end>=rcEnd){
				leftOver->push_back(rcEnd);
			}
			if(dogMode){
				kmer seq(str2num(line.substr(0,anchorSize))),rcSeq(rcb(seq,anchorSize)),canon(min(seq,rcSeq));
				anchors->push_back(canon);
				for(uint j(0);j+anchorSize<line.size();++j){
					updateK(seq,line[j+anchorSize]);
					updateRCK(rcSeq,line[j+anchorSize]);
					if((j+1)%fracKmer==0){
						canon=(min(seq, rcSeq));
						anchors->push_back(canon);
					}
				}
			}
		}
	}
	auto end1=system_clock::now();auto waitedFor1=end1-start1;cout<<"Duration "<<duration_cast<seconds>(waitedFor1).count()<<" seconds"<<endl;

	cout<<"Sorting anchors: "<<flush;auto start2=system_clock::now();
	sort( leftOver->begin(), leftOver->end() );
	leftOver->erase( unique( leftOver->begin(), leftOver->end() ), leftOver->end() );
	sort( rightOver->begin(), rightOver->end() );
	rightOver->erase( unique( rightOver->begin(), rightOver->end() ), rightOver->end() );
	auto end2=system_clock::now();auto waitedFor2=end2-start2;cout<<"Duration "<<duration_cast<seconds>(waitedFor2).count()<<" seconds"<<endl;

	cout<<"Creating MPHF: "<<flush;auto start3=system_clock::now();
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
	auto end3=system_clock::now();auto waitedFor3=end3-start3;cout<<"Duration "<<duration_cast<seconds>(waitedFor3).count()<<" seconds"<<endl;

	cout<<"Filling index: "<<flush;auto start4=system_clock::now();
	leftIndicesstr.resize(leftsize,{});
	rightIndicesstr.resize(rightsize,{});
	anchorsChecking.resize(anchorNumber,0);
	fillIndicesstrbutanchors();
	auto end4=system_clock::now();auto waitedFor4=end4-start4;cout<<"Duration "<<duration_cast<seconds>(waitedFor4).count()<<" seconds"<<endl;
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
					if(anchorsPositionVector[hash].size()<maxPositionAnchors){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({i,0});
						mutexV[hash%1000].unlock();
					}
				}else{
					anchorsPosition[hash]={i,0};
				}
			}else{
				if(vectorMode){
					if(anchorsPositionVector[hash].size()<maxPositionAnchors){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({-i,0});
						mutexV[hash%1000].unlock();
					}
				}else{
					anchorsPosition[hash]={-i,0};
				}
			}
			anchorsChecking[hash]=getHash(canon);
			for(uint j(0);j+anchorSize<line.size();++j){
				updateK(seq,line[j+anchorSize]);
				updateRCK(rcSeq,line[j+anchorSize]);
				if((j+1)%fracKmer==0){
					canon=(min(seq, rcSeq));
					uint64_t hash=anchorsMPHF.lookup(canon);
					if(canon==seq){
						if(vectorMode){
							if(anchorsPositionVector[hash].size()<maxPositionAnchors){
								mutexV[hash%1000].lock();
								anchorsPositionVector[hash].push_back({i,j+1});
								mutexV[hash%1000].unlock();
							}
						}else{
							anchorsPosition[hash]={i,j+1};
						}
					}else{
						if(vectorMode){
							if(anchorsPositionVector[hash].size()<maxPositionAnchors){
								mutexV[hash%1000].lock();
								anchorsPositionVector[hash].push_back({-i,j+1});
								mutexV[hash%1000].unlock();
							}
						}else{
							anchorsPosition[hash]={-i,j+1};
						}
					}
					anchorsChecking[hash]=getHash(canon);
				}
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
		}
		if(beg>=rcBeg){
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
		}
		if(end>=rcEnd){
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
					if(anchorsPositionVector[hash].size()<maxPositionAnchors){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({i,0});
						mutexV[hash%1000].unlock();
					}
				}else{
					anchorsPosition[hash]={i,0};
				}
			}else{
				if(vectorMode){
					if(anchorsPositionVector[hash].size()<maxPositionAnchors){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({-i,0});
						mutexV[hash%1000].unlock();
					}
				}else{
					anchorsPosition[hash]={-i,0};
				}
			}
			anchorsChecking[hash]=getHash(canon);
			for(uint j(0);j+anchorSize<line.size();++j){
				seq=seq.substr(1)+line[j+anchorSize];
				rcSeq=revCompChar(line[j+anchorSize])+rcSeq.substr(0,anchorSize-1);
				if((j+1)%fracKmer==0){
					canon=(min(seq,rcSeq));
					int64_t hash=anchorsMPHFstr.lookup(canon);
					if(canon==seq){
						if(vectorMode){
							if(anchorsPositionVector[hash].size()<maxPositionAnchors){
								mutexV[hash%1000].lock();
								anchorsPositionVector[hash].push_back({i,j+1});
								mutexV[hash%1000].unlock();
							}
						}else{
							anchorsPosition[hash]={i,j+1};
						}
					}else{
						if(vectorMode){
							if(anchorsPositionVector[hash].size()<maxPositionAnchors){
								mutexV[hash%1000].lock();
								anchorsPositionVector[hash].push_back({-i,j+1});
								mutexV[hash%1000].unlock();
							}
						}else{
							anchorsPosition[hash]={-i,j+1};
						}
					}
					anchorsChecking[hash]=getHash(canon);
				}
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
		}
		if(beg>=rcBeg){
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
		}
		if(end>=rcEnd){
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
					if(anchorsPositionVector[hash].size()<maxPositionAnchors){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({i,0});
						mutexV[hash%1000].unlock();
					}
				}else{
					anchorsPosition[hash]={i,0};
				}
			}else{
				if(vectorMode){
					if(anchorsPositionVector[hash].size()<maxPositionAnchors){
						mutexV[hash%1000].lock();
						anchorsPositionVector[hash].push_back({-i,0});
						mutexV[hash%1000].unlock();
					}
				}else{
					anchorsPosition[hash]={-i,0};
				}
			}
			anchorsChecking[hash]=getHash(canon);
			for(uint j(0);j+anchorSize<line.size();++j){
				updateK(seq,line[j+anchorSize]);
				updateRCK(rcSeq,line[j+anchorSize]);
				if((j+1)%fracKmer==0){
					canon=(min(seq, rcSeq));
					uint64_t hash=anchorsMPHF.lookup(canon);
					if(canon==seq){
						if(vectorMode){
							if(anchorsPositionVector[hash].size()<maxPositionAnchors){
								mutexV[hash%1000].lock();
								anchorsPositionVector[hash].push_back({i,j+1});
								mutexV[hash%1000].unlock();
							}
						}else{
							anchorsPosition[hash]={i,j+1};
						}
					}else{
						if(vectorMode){
							if(anchorsPositionVector[hash].size()<maxPositionAnchors){
								mutexV[hash%1000].lock();
								anchorsPositionVector[hash].push_back({-i,j+1});
								mutexV[hash%1000].unlock();
							}
						}else{
							anchorsPosition[hash]={-i,j+1};
						}
					}
					anchorsChecking[hash]=getHash(canon);
				}
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
		}
		if(beg>=rcBeg){
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
		}
		if(end>=rcEnd){
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
	if(not unitigFile.good()){
		cout<<"Problem with unitigs file"<<endl;
		exit(0);
	}
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
	sizeBuckets=unitigs.size()/(nbBuckets-1);
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
			delete(readFile);
			//~ readFile.open(file);
			readFile=new zstr::ifstream(file);
			//~ readFileF=fopen(file.c_str(),"r");
			cout<<"Mapping: "<<file<<" on "<<unitigFileName<<endl;
			if(not readFile->good()){
				cout<<"A probleme with file: "<<file<<endl;
			}
			unsigned char nbThreads(coreNumber);
			vector<thread> threads;
			for (size_t i(0); i<nbThreads; ++i){
				threads.push_back(thread(&Aligner::alignPartGreedy,this,i));
			}
			for(auto &t : threads){t.join();}
			last=i+1;
		}
	}
	file=reads.substr(last);
	delete(readFile);
	cout<<file<<endl;
	readFile=new zstr::ifstream(file);
	//~ readFile.open(file);
	if(not readFile->good()){
		cout<<"A probleme with file: "<<file<<endl;
		return;
	}
	//~ readFileF=fopen(file.c_str(),"r");
	cout<<"Mapping: "<<file<<" on "<<unitigFileName<<endl;
	unsigned char nbThreads(coreNumber);
	vector<thread> threads;
	for (size_t i(0); i<nbThreads; ++i){
		threads.push_back(thread(&Aligner::alignPartGreedy,this,i));
	}
	for(auto &t : threads){t.join();}

	cout<<"The End"<<endl;
	cout<<"Reads: "<<intToString(readNumber)<<endl;
	cout<<"Not anchored : "<<intToString(noOverlapRead)<<" Percent: "<<(100*float(noOverlapRead))/readNumber<<endl;
	cout<<"Anchored and aligned : "<<intToString(alignedRead)<<" Percent: "<<(100*float(alignedRead))/(alignedRead+notAligned)<<endl;
	cout<<"Not aligned: "<<intToString(notAligned)<<" Percent: "<<(100*float(notAligned))/(alignedRead+notAligned)<<endl;
	auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
	cout<<"Reads/seconds: "<<intToString(readNumber/(chrono::duration_cast<chrono::seconds>(waitedFor).count()+1))<<endl;
	cout<<"Mapping in seconds: "<<intToString(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<endl;
	if(pairedMode){
		cout<<"Super reads: "<<intToString(superReads)<<endl;
		cout<<"Failed super reads compaction: "<<intToString(notCompatedSR)<<endl;
		cout<<"Percentage paired: "<<(float(100)*superReads/(superReads+notCompatedSR))<<endl;
		cout<<"Included path: "<<intToString(includedPath)<<" Overlapping path: "<<intToString(overlappingPath)<<" Overlapping Str: "<<intToString(overlappingStr)<<"  Single unitig "<<intToString(singleMiddle)<<" Neighbors "<<intToString(neighbor)<<endl;
		cout<<"Meaningful pair-mapping rate: "<<overlappingPath+overlappingStr+singleMiddle+neighbor<<"wrt"<<failed_pair<<"="<<100*(overlappingPath+overlappingStr+singleMiddle+neighbor)/(overlappingPath+overlappingStr+singleMiddle+neighbor+failed_pair)<<"%"<<endl;

	}
	cout<<endl;
}



string Aligner::printPath(const vector<int32_t>& path){
	string res;
	for(uint i(0); i<path.size(); ++i){
		res+=to_string(path[i]);
		res+='.';
	}
	res+='\n';
	return res;
}


string Aligner::path2nuc(const vector<int32_t>& path){
	string res;
	res+=to_string(path[1])+":";//ANCHOR
	for(int i(2);i<path.size();++i){
		if(path[i]>0){
			res+=unitigs[path[i]][0];
		}else{
			res+=(((unitigsRC[-path[i]]))[0]);
		}
	}
	res+=":";
	return res;
}
















