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



#include <thread>
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



uint trimingBases(0);



//TODO, SEARCH FOR FEW ANCHORS, TRY TO MAP, IF NOT MAPPED SEARCH FOR MORE
vector<uNumber> Aligner::alignReadGreedy(const string& read, bool& overlapFound, uint errors, bool& rc){
	vector<pair<kmer,uint>> listOverlap(getNOverlap(read,tryNumber));
	if(listOverlap.empty()){++noOverlapRead;return {};}
	overlapFound=true;
	vector<uNumber> pathBegin,pathEnd;
	for(uint start(0); start<(uint)listOverlap.size(); ++start){
		pathBegin={};
		uint errorBegin(checkBeginGreedy(read,listOverlap[start],pathBegin,errors));
		if(errorBegin<=errors){
			pathEnd={};
			uint errorsEnd(checkEndGreedy(read,listOverlap[start],pathEnd,errors-errorBegin));
			if(errorsEnd+errorBegin<=errors){
				//~ ++alignedRead;
				reverse(pathBegin.begin(),pathBegin.end());
				pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
				return pathBegin;
			}
		}
	}
	//~ if(!rc){rc=true;return alignReadGreedy(reverseComplements(read), overlapFound,errors, rc);}
	//~ ++notAligned;
	return {};
}



vector<uNumber> Aligner::alignReadGreedyAnchors(const string& read, uint errorMax,const pair<pair<uint,uint>,uint>& anchor,uint& errors){
	vector<uNumber> pathBegin,pathEnd;
	string unitig("");
	bool returned(false);
	int unitigNumber(anchor.first.first),positionUnitig(anchor.first.second),positionRead(anchor.second);
	if(unitigNumber>=0){
		unitig=(unitigs[unitigNumber]);
	}else{
		if(rcMode){
			unitig=(unitigsRC[-unitigNumber]);
		}else{
			unitig=reverseComplements(unitigs[-unitigNumber]);
		}
		positionUnitig=unitig.size()-positionUnitig-anchorSize;
		returned=true;
	}
	if(positionRead>=positionUnitig){
		if(read.size()-positionRead>=unitig.size()-positionUnitig){
			//CASE 1 : unitig included in read
			//~ cout<<"1:"<<endl;
			errors=(missmatchNumber(read.substr(positionRead-positionUnitig,unitig.size()),unitig,errorMax));
			if(errors<=errorMax){
				pathBegin={};
				//~ uint errorBegin;
				errors+=(checkBeginGreedy(read,{str2num(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
				if(errors<=errorMax){
					pathEnd={(int)unitigNumber};
					//~ uint errorsEnd;
					errors+=(checkEndGreedy(read,{str2num(unitig.substr(unitig.size()-k+1,k-1)),positionRead-positionUnitig+unitig.size()},pathEnd,errorMax-errors));
					if(errors<=errorMax){
						reverse(pathBegin.begin(),pathBegin.end());
						pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
						return pathBegin;
					}
				}
			}
		}else{
			//CASE 2 : unitig overap read
			//~ cout<<"2:"<<endl;
			errors=(missmatchNumber(read.substr(positionRead-positionUnitig),unitig.substr(0,read.size()-positionRead+positionUnitig),errorMax));
			if(errors<=errorMax){
				pathBegin={};
				//~ uint errorBegin;
				errors+=(checkBeginGreedy(read,{str2num(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
				if(errors<=errorMax){
					reverse(pathBegin.begin(),pathBegin.end());
					pathBegin.push_back((int)unitigNumber);
					return pathBegin;
				}
			}
		}
	}else{
		if(read.size()-positionRead>=unitig.size()-positionUnitig){
			//CASE 3 : read overlap unitig
			//~ cout<<"3:"<<endl;
			errors=(missmatchNumber(unitig.substr(positionUnitig-positionRead),read.substr(0,unitig.size()+positionRead-positionUnitig),errorMax));
			if(errors<=errorMax){
				pathEnd={(int)positionUnitig-(int)positionRead,(int)unitigNumber};
				//~ uint errorsEnd;
				errors+=(checkEndGreedy(read,{str2num(unitig.substr(unitig.size()-k+1,k-1)),positionRead-positionUnitig+unitig.size()},pathEnd,errorMax-errors));
				if(errors<=errorMax){
					return pathEnd;
				}
			}
		}else{
			//CASE 4 : read included in unitig
			//~ cout<<"4:"<<endl;
			//~ cout<<unitig.substr(positionUnitig-positionRead,read.size())<<endl;
			errors=(missmatchNumber(unitig.substr(positionUnitig-positionRead,read.size()),read,errorMax));
			if(errors<=errorMax){
				return {(int)positionUnitig-(int)positionRead,(int)unitigNumber};
			}else{
				//~ cout<<read<<endl;
				//~ cin.get();
			}
		}
	}
	return {};
}


//TODO SEGFAULT HERE
vector<uNumber> Aligner::alignReadGreedyAnchorsstr(const string& read, uint errorMax, const pair<pair<uint,uint>,uint>& anchor, uint& errors){
	//~ cout<<"go"<<endl;
	vector<uNumber> pathBegin,pathEnd;
	string unitig("");
	bool returned(false);
	int unitigNumber(anchor.first.first),positionUnitig(anchor.first.second),positionRead(anchor.second);
	if(unitigNumber>=0){
		unitig=(unitigs[unitigNumber]);
	}else{
		if(rcMode){
			unitig=(unitigsRC[-unitigNumber]);
		}else{
			unitig=reverseComplements(unitigs[-unitigNumber]);
		}
		positionUnitig=unitig.size()-positionUnitig-anchorSize;
		returned=true;
	}

	if(positionRead>=positionUnitig){
		if(read.size()-positionRead>=unitig.size()-positionUnitig){
			//~ cout<<1<<endl;
			errors=(missmatchNumber(read.substr(positionRead-positionUnitig,unitig.size()),unitig,errorMax));
			if(errors<=errorMax){
				pathBegin={};
				errors+=(checkBeginGreedy(read,{(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
				if(errors<=errorMax){
					pathEnd={(int)unitigNumber};
					errors+=(checkEndGreedy(read,{(unitig.substr(unitig.size()-k+1,k-1)),positionRead-positionUnitig+unitig.size()},pathEnd,errorMax-errors));
					if(errors<=errorMax){
						reverse(pathBegin.begin(),pathBegin.end());
						pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
						return pathBegin;
					}
				}
			}
		}else{
			//~ cout<<2<<endl;
			errors=(missmatchNumber(read.substr(positionRead-positionUnitig),unitig.substr(0,read.size()-positionRead+positionUnitig),errorMax));
			if(errors<=errorMax){
				pathBegin={};
				errors+=(checkBeginGreedy(read,{(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
				if(errors<=errorMax){
					reverse(pathBegin.begin(),pathBegin.end());
					pathBegin.push_back((int)unitigNumber);
					return pathBegin;
				}
			}
		}
	}else{
		if(read.size()-positionRead>=unitig.size()-positionUnitig){
			//~ cout<<3<<endl;
			errors=(missmatchNumber(unitig.substr(positionUnitig-positionRead),read.substr(0,unitig.size()+positionRead-positionUnitig),errorMax));
			if(errors<=errorMax){
				pathEnd={(int)positionUnitig-(int)positionRead,(int)unitigNumber};
				errors+=(checkEndGreedy(read,{unitig.substr(unitig.size()-k+1,k-1),positionRead-positionUnitig+unitig.size()},pathEnd,errorMax-errors));
				if(errors<=errorMax){
					return pathEnd;
				}
			}
		}else{
			//~ cout<<4<<endl;
			uint errors(missmatchNumber(unitig.substr(positionUnitig-positionRead,read.size()),read,errorMax));
			if(errors<=errorMax){
				return {(int)positionUnitig-(int)positionRead,(int)unitigNumber};
			}
		}
	}
	return {};
}



uint Aligner::mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	//~ cout<<"moleg"<<endl;
	if(overlap.second<=trimingBases){path.push_back(0);return 0;}
	string unitig,readLeft(read.substr(0,overlap.second)),nextUnitig;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getEnd(overlap.first));
	uint miniMiss(1000),miniMissIndice(9);
	bool ended(false),equal(false);
	int offset(0);
	kmer nextOverlap(0);
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				path.push_back(unitig.size()-readLeft.size()-k+1);
				return 0;
			}else if(miss==miniMiss){
				equal=true;
			}else if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
				equal=false;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnLeftEndGreedy(read , path, {str2num(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)}, errors);

			}else if(miss==miniMiss){
				equal=true;
			}else if(miss<miniMiss){
				kmer overlapNum(str2num(unitig.substr(0,k-1)));
				if(miss<miniMiss){
					ended=false;
					equal=false;
					miniMiss=miss;
					miniMissIndice=i;
					nextUnitig=unitig;
					nextOverlap=overlapNum;
				}
			}
		}
	}

	if (miniMiss<=errors and not equal){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(ended){
			path.push_back(offset);
			return miniMiss;
		}
		return miniMiss+mapOnLeftEndGreedy(read , path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)}, errors-miniMiss);
	}
	if(equal){
		return 1000;
	}
	return miniMiss;
}



uint Aligner::mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<string, uint>& overlap , uint errors){
	if(overlap.second<=trimingBases){path.push_back(0);return 0;}
	string unitig,readLeft(read.substr(0,overlap.second)),nextUnitig;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getEnd(overlap.first));
	uint miniMiss(1000),miniMissIndice(9);
	bool ended(false),equal(false);
	int offset(0);
	string nextOverlap;
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				path.push_back(unitig.size()-readLeft.size()-k+1);
				return 0;
			}else if(miss==miniMiss){
				equal=true;
			}else if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
				equal=false;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnLeftEndGreedy(read , path, {(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)}, errors);
			}else if(miss==miniMiss){
				equal=true;
			}else if(miss<miniMiss){
				ended=false;
				equal=false;
				miniMiss=miss;
				miniMissIndice=i;
				nextUnitig=unitig;
				nextOverlap=((unitig.substr(0,k-1)));
			}
		}
	}
	if (miniMiss<=errors and not equal){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(ended){
			path.push_back(offset);
			return miniMiss;
		}
		return miniMiss+mapOnLeftEndGreedy(read , path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)}, errors-miniMiss);
	}
	if(equal){
		return 1000;
	}
	return miniMiss;
}



uint Aligner::mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	//~ cout<<"moreg"<<endl;
	string unitig,readLeft(read.substr(overlap.second)),nextUnitig;
	if(readLeft.size()<trimingBases){return 0;}
	//~ cout<<readLeft<<endl;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=getBegin(overlap.first);
	uint miniMiss(1001), miniMissIndice(9);
	bool ended(false),equal(false);
	kmer nextOverlap(0);
	//~ cout<<rangeUnitigs.size()<<endl;
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);

		//case the rest of the read is too small
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			//~ cout<<unitig.substr(k-1,readLeft.size())<<endl;
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return 0;
			}else if(miss==miniMiss){
				equal=true;
			}else if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
				equal=false;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
			//~ cout<<unitig.substr(k-1)<<endl;
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return (mapOnRightEndGreedy(read , path, {str2num(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)}, errors));
			}else if(miss==miniMiss){
				equal=true;
			}else if(miss<miniMiss){
				ended=false;
				equal=false;
				kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
				miniMiss=miss;
				miniMissIndice=i;
				nextUnitig=unitig;
				nextOverlap=overlapNum;
			}
		}
	}
	if(miniMiss<=errors and not equal){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(ended){return miniMiss;}
		return (miniMiss+mapOnRightEndGreedy(read , path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)}, errors-miniMiss));
	}
	if(equal){
		return 1000;
	}
	return miniMiss;
}



uint Aligner::mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<string, uint>& overlap , uint errors){
	string unitig,readLeft(read.substr(overlap.second)),nextUnitig;
	if(readLeft.size()<trimingBases){return 0;}
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=getBegin(overlap.first);
	uint miniMiss(1000), miniMissIndice(9);
	bool ended(false),equal(false);
	string nextOverlap;
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return 0;
			}else if(miss==miniMiss){
				equal=true;
			}else if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
				equal=false;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return (mapOnRightEndGreedy(read , path, {(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)}, errors));
			}else if(miss==miniMiss){
				equal=true;
			}else if(miss<miniMiss){
				ended=false;
				equal=false;
				miniMiss=miss;
				miniMissIndice=i;
				nextUnitig=unitig;
				nextOverlap=(unitig.substr(unitig.size()-k+1,k-1));
			}
		}
	}
	if(miniMiss<=errors and not equal){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(ended){return miniMiss;}
		return (miniMiss+mapOnRightEndGreedy(read , path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)}, errors-miniMiss));
	}
	if(equal){
		return 1000;
	}
	return miniMiss;
}



uint Aligner::checkBeginGreedy(const string& read,const pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors){
	//~ cout<<"cbg"<<endl;
	if(overlap.second<=trimingBases){path.push_back(0);return 0;}
	string readLeft(read.substr(0,overlap.second)),unitig,nextUnitig;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getEnd(overlap.first));
	uint minMiss(1001),indiceMinMiss(9);
	bool ended(false),equal(false);
	int offset(0);
	kmer nextOverlap(0);
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=rangeUnitigs[i].first;
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				path.push_back(unitig.size()-readLeft.size()-k+1);
				return 0;
			}else if(miss==minMiss){
				equal=true;
			}else if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				equal=false;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnLeftEndGreedy(read, path, {str2num(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)},errors);
			}else if(miss==minMiss){
				equal=true;
			}else if(miss<minMiss){
				kmer overlapNum(str2num(unitig.substr(0,k-1)));
				ended=false;
				equal=false;
				minMiss=miss;
				indiceMinMiss=i;
				nextUnitig=unitig;
				nextOverlap=overlapNum;
			}

		}
	}
	if(minMiss<=errors and not equal){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){
			path.push_back(offset);
			return minMiss;
		}
		return minMiss+mapOnLeftEndGreedy(read, path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)},errors-minMiss);
	}
	if(equal){
		return 1000;
	}
	return minMiss;
}



uint Aligner::checkBeginGreedy(const string& read,const pair<string, uint>& overlap, vector<uNumber>& path, uint errors){
	if(overlap.second<=trimingBases){path.push_back(0);return 0;}
	string readLeft(read.substr(0,overlap.second)),unitig,nextUnitig;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getEnd(overlap.first));
	uint minMiss(1001),indiceMinMiss(9);
	bool ended(false),equal(false);
	int offset(0);
	string nextOverlap;
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=rangeUnitigs[i].first;
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				path.push_back(unitig.size()-readLeft.size()-k+1);
				return 0;
			}else if(miss==minMiss){
				equal=true;
			}else if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				equal=false;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnLeftEndGreedy(read, path, {(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)},errors);
			}else if(miss==minMiss){
				equal=true;
			}else if(miss<minMiss){
				ended=false;
				equal=false;
				minMiss=miss;
				indiceMinMiss=i;
				nextUnitig=unitig;
				nextOverlap=((unitig.substr(0,k-1)));
			}

		}
	}
	if(minMiss<=errors and not equal){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){
			path.push_back(offset);
			return minMiss;
		}
		return minMiss+mapOnLeftEndGreedy(read, path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)},errors-minMiss);
	}
	if(equal){
		return 1000;
	}
	return minMiss;
}



uint Aligner::checkEndGreedy(const string& read, const pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors){
	//~ cout<<"ceg"<<endl;
	string readLeft(read.substr(overlap.second)),unitig,nextUnitig;
	if(readLeft.size()<=trimingBases){return 0;}
	//~ cout<<readLeft<<endl;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getBegin(overlap.first));
	uint minMiss(1001),indiceMinMiss(9);
	bool ended(false),equal(false);
	kmer nextOverlap(0);
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=rangeUnitigs[i].first;
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			//~ cout<<unitig<<endl;
			//~ cout<<"end"<<endl;
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return 0;
			}else if(miss==minMiss){
				equal=true;
			}else if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				equal=false;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
			//~ cout<<unitig<<endl;
			//~ cout<<"noend"<<endl;
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnRightEndGreedy(read, path, {str2num(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)},errors);
			}else if(miss==minMiss){
				equal=true;
			}else if(miss<minMiss){
				kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
				minMiss=miss;
				indiceMinMiss=i;
				nextOverlap=overlapNum;
				nextUnitig=unitig;
				ended=false;
				equal=false;
			}
		}
	}
	if(minMiss<=errors and not equal){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){return minMiss;}
		return minMiss+mapOnRightEndGreedy(read, path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)},errors-minMiss);
	}
	if(equal){
		return 1000;
	}
	return minMiss;
}



uint Aligner::checkEndGreedy(const string& read, const pair<string, uint>& overlap, vector<uNumber>& path, uint errors){
	string readLeft(read.substr(overlap.second)),unitig,nextUnitig;
	if(readLeft.size()<=trimingBases){return 0;}
	vector<pair<string,uNumber>> rangeUnitigs(getBegin(overlap.first));
	uint minMiss(1000),indiceMinMiss(9);
	bool ended(false),equal(false);
	string nextOverlap;
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=rangeUnitigs[i].first;
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return 0;
			}
			else if(miss==minMiss){
				equal=true;
			}else if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				equal=false;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnRightEndGreedy(read, path, {(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)},errors);
			}else if(miss==minMiss){
				equal=true;
			}else if(miss<minMiss){
				nextOverlap=((unitig.substr(unitig.size()-k+1,k-1)));
				minMiss=miss;
				indiceMinMiss=i;
				nextUnitig=unitig;
				equal=ended=false;
			}
		}
	}
	if(minMiss<=errors and not equal){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){return minMiss;}
		return minMiss+mapOnRightEndGreedy(read, path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)},errors-minMiss);
	}
	if(equal){
		return 1000;
	}
	return minMiss;
}



vector<int> Aligner::inclued(vector<int>& v1, vector<int>& v2){
	if(v1.size()<v2.size()){
		for(uint i(1);i<v1.size();++i){
			bool found (false);
			for(uint j(1);j<v2.size();++j){
				if(v1[i]==v2[j] or v1[i]==-v2[j]){
					found=true;
					break;
				}
			}
			if(not found){
				//~ for(uint iii(1);iii<v1.size();++iii){
					//~ cout<<v1[iii]<<" ";
				//~ }
				//~ cout<<endl;
				//~ for(uint iii(1);iii<v2.size();++iii){
					//~ cout<<v2[iii]<<" ";
				//~ }
				//~ cout<<endl;
				//~ cin.get();
				//~ return {};
			}
		}
		return v1;
	}else{
		for(uint i(1);i<v2.size();++i){
			bool found (false);
			for(uint j(1);j<v1.size();++j){
				if(v2[i]==v1[j] or v2[i]==-v1[j]){
					found=true;
					break;
				}
			}
			if(not found){
				//~ for(uint iii(1);iii<v1.size();++iii){
					//~ cout<<v1[iii]<<" ";
				//~ }
				//~ cout<<endl;
				//~ for(uint iii(1);iii<v2.size();++iii){
					//~ cout<<v2[iii]<<" ";
				//~ }
				//~ cout<<endl;
				//~ cin.get();
				return {};
			}
		}
		return v2;
	}

	return {};
}



void Aligner::alignReadOpti(const string& read, vector<int>& path,bool perfect=false){
	path={};
	vector<int> pathMem;
	uint errors(0);
	vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	vector<int> errorsFromPreviousMapping(listAnchors.size(),-1);
	if(listAnchors.empty()){
		++noOverlapRead;
		return;
	}
	string superRead,superReadMem;
	//~ random_shuffle ( listAnchors.begin(), listAnchors.end() );
	while(errors<=errorsMax){
		bool found(false);
		pathMem={};
		uint errorInMapping(0);
		for(uint i(0);i<listAnchors.size();++i){
			path={};
			errorInMapping=0;
			if(errorsFromPreviousMapping[i]<=(int)errors){
				if(stringMode){
					path=alignReadGreedyAnchorsstr(read,errors,listAnchors[i],errorInMapping);
				}else{
					path=alignReadGreedyAnchors(read,errors,listAnchors[i],errorInMapping);
				}
				errorsFromPreviousMapping[i]=errorInMapping;
				//MAPPING IS FOUND
				if(not path.empty()){
					if(noMultiMapping){
						if(found){
							superRead=(recoverSuperReadsCor(path,read.size()));
							superReadMem=(recoverSuperReadsCor(pathMem,read.size()));
							if(superRead!=superReadMem){
								//~ cout<<superRead<<endl;
								//~ cout<<superReadMem<<endl;
								//~ cin.get();
								path={};
								return;
							}else{
								if(path.size()<pathMem.size()){
									pathMem=path;
								}
							}
						}else{
							pathMem=path;
						}
						found=true;
					}else{
						//DONE
						++alignedRead;
						return;
					}
				}
			}else{
			}
		}
		++errors;
		if(perfect or found){
			path=pathMem;
			++alignedRead;
			return;
		}
	}
	path=pathMem;
	if(not path.empty()){
		++alignedRead;
	}else{
		//~ cout<<">1"<<"\n"<<read<<endl;
	}
}




void Aligner::alignReadAllOpti(const string& read, vector<vector<int>>& pathVector){
	pathVector={};
	vector<int> path;
	uint errors(0);
	vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	if(listAnchors.empty()){++noOverlapRead;return;}
	while(errors<=errorsMax){
		bool found(false);
		uint errorInMapping(0);
		for(uint i(0);i<listAnchors.size();++i){
			path={};
			if(stringMode){
				path=alignReadGreedyAnchorsstr(read,errors,listAnchors[i],errorInMapping);
			}else{
				path=alignReadGreedyAnchors(read,errors,listAnchors[i],errorInMapping);
			}
			//MAPPING IS FOUND
			if(not path.empty()){
				pathVector.push_back(path);
				found=true;
			}
		}
		++errors;
		if(found){
			++alignedRead;
			return;
		}
	}
}



void Aligner::alignReadAll(const string& read, vector<vector<int>>& pathVector){
	pathVector={};
	vector<int> path;
	uint errors(0);
	vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	if(listAnchors.empty()){++noOverlapRead;return;}
	while(errors<=errorsMax){
		for(uint i(0);i<listAnchors.size();++i){
			path={};
			uint errorInMapping(0);
			if(stringMode){
				path=alignReadGreedyAnchorsstr(read,errors,listAnchors[i],errorInMapping);
			}else{
				path=alignReadGreedyAnchors(read,errors,listAnchors[i],errorInMapping);
			}
			//MAPPING IS FOUND
			if(not path.empty()){
				pathVector.push_back(path);
			}
		}
		++errors;
	}
}



void Aligner::alignReadFrom(const string& read, vector<int>& path, int unumber){
	vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	for(uint i(0);i<listAnchors.size();++i){
		auto  anchor=listAnchors[i];
		//~ if((anchor.first.first==unumber or anchor.first.first==-unumber) and anchor.second==0 and anchor.first.second==0){
			path={};
			uint errorInMapping(0);
			if(stringMode){
				path=alignReadGreedyAnchorsstr(read,0,anchor,errorInMapping);
			}else{
				path=alignReadGreedyAnchors(read,0,anchor,errorInMapping);
			}
			if(not path.empty()){
				return;
			}

		//~ }
	}
}




void Aligner::alignPartGreedy(uint indiceThread){
	vector<pair<string,string>> multiread;
	vector<uNumber> path,path2;
	vector<vector<int>> pathVector;
	string read,read2,header,header2,corrected,superRead,toWrite;
	vector<string> toWriteComp(nbBuckets);
	pair<string,string> superpath;
	while(not readFile.eof()){
	//~ while(not feof(readFileF)){
		toWrite="";
		if(keepOrder){
			while(threadToRead!=indiceThread){
				this_thread::sleep_for (chrono::microseconds(1));
			}
			readMutex.lock();
			getReads(multiread,100000);
			threadToRead=(threadToRead+1)%coreNumber;
			readMutex.unlock();
		}else{
			readMutex.lock();
			getReads(multiread,100000);
			readMutex.unlock();
		}
		if(pairedMode){
			//PAIRED READS
			for(uint i(0);i+1<multiread.size();i+=2){
				header=multiread[i].first;
				read=multiread[i].second;
				header2=multiread[i+1].first;
				read2=multiread[i+1].second;
				readNumber+=2;
				bool rc(false), noOverlap(false), overlapFound(false);
				alignReadOpti(read,path);
				alignReadOpti(read2,path2);
				if(path.empty()){
					++notAligned;
				}else{
					//~ path=cleanSR(path,read.size());
				}
				if(path2.empty()){
					++notAligned;
				}else{
					//~ path2=cleanSR(path2,read2.size());
				}
				superpath=(recoverSuperReadsPairedNoStr(path,path2));
				if(superpath.first!=""){
					if(superpath.second!=""){
						if(headerNeeded){
							toWrite+=header+'\n'+superpath.first+'\n'+header2+'\n'+superpath.second+'\n';
						}else{
							toWrite+=superpath.first+'\n'+superpath.second+'\n';
						}
					}else{
						if(headerNeeded){
							toWrite+=header+'\n'+superpath.first+'\n';
						}else{
							toWrite+=superpath.first+'\n';
						}
					}
				}else{
					if(superpath.second!=""){
						if(headerNeeded){
							toWrite+=header2+'\n'+superpath.second+'\n';
						}else{
							toWrite+=superpath.second+'\n';
						}
					}
				}
			}
		}else{
			//UNPAIRED READS
			for(uint i(0);i<multiread.size();i++){
				header=multiread[i].first;
				read=multiread[i].second;
				++readNumber;
				if(uniqueOptimalMappingMode){
					alignReadOpti(read,path);
					if(path.empty()){
						if(correctionMode){
							toWrite+=header+'\n'+read+'\n';
						}
						++notAligned;
						continue;
					}
					if(correctionMode){
						//CORRECTION
						if(not path.empty()){
							superRead=(recoverSuperReadsCor(path,read.size()));
							if(superRead!=""){
								toWrite+=header+'\n'+superRead+'\n';
							}else{
								toWrite+=header+'\n'+read+'\n';
							}
						}else{
							//~ toWrite+=header+'\n'+read+'\n';
						}
					}else if(preciseOutput){
						//PRECISE MODE
						uint position(path[0]);
						path=vector<uNumber>(&path[1],&path[path.size()]);
						superRead=(recoverSuperReads(path));
						superRead=superRead.substr(position);
						path.push_back(position);
						int lastUnitigNumber(path[path.size()-2]);
						if (lastUnitigNumber<0){
							lastUnitigNumber=-lastUnitigNumber;
						}
						path.push_back((int)unitigs[lastUnitigNumber].size()+(int)read.size()-(int)superRead.size());
						superRead=(recoverSuperReadsNoStr(path,0));
						if(superRead!=""){
							toWrite+=header+'\n'+superRead+'\n';
						}
					}else if(compressionMode){
						//~ cout<<"test"<<endl;
						//~ cout<<path[1]<<endl;
						//~ cout<<sizeBuckets<<endl;
						int goodBucket(path[1]);
						if(goodBucket<0){goodBucket=-goodBucket;}
						goodBucket/=sizeBuckets;
						//~ cout<<goodBucket<<endl;

						//~ cout<<goodBucket<<endl;
						toWriteComp[goodBucket]+=path2nuc(path);//WRITE ANCHORS AND NUCLEOTIDE TO MAKE THE PATH
						//~ toWriteComp[goodBucket]+=to_string(path[0])+":"+to_string(read.size())+":";//Read position and size
						toWriteComp[goodBucket]+=to_string(path[0])+":";//Read position
						superRead=(recoverSuperReadsCor(path,read.size()));
						//~ toWriteComp[goodBucket]+=codeMiss(read, superRead);//ENCODE THE MISSMATCHES
						toWriteComp[goodBucket]+=(char)255;
						//~ cout<<"end"<<endl;
					}else{
						//REGULAR MODE
						path=cleanSR(path,read.size());
						superRead=(recoverSuperReadsNoStr(path,1));
						if(superRead!=""){
							if(printAlignment){
								toWrite+=header+'\n'+superRead+'\n';
								toWrite+=read+'\n'+recoverSuperReadsCor(path,read.size())+'\n';
							}else{
								if(headerNeeded){
									toWrite+=superRead+'\n';
								}else{
									toWrite+=header+'\n'+superRead+'\n';
								}
							}
						}
					}
				}else{
					if(optimalMappingMode){
						alignReadAllOpti(read,pathVector);
					}else{
						alignReadAll(read,pathVector);
					}
					if(pathVector.empty()){++notAligned;}else{++alignedRead;}
					sort(pathVector.begin(), pathVector.end());
					pathVector.erase(unique(pathVector.begin(), pathVector.end()), pathVector.end());
					for(uint i(0);i<pathVector.size();++i){
						superRead=(recoverSuperReadsNoStr(pathVector[i],1));
						if(superRead!=""){
							toWrite+=header+'\n'+superRead+'\n';
							if(printAlignment){
								toWrite+=read+'\n'+recoverSuperReadsCor(pathVector[i],read.size())+'\n';
							}
						}else{
							cout<<"wut"<<endl;
						}
					}
				}
			}
		}
		if(keepOrder){
			while(threadToPrint!=indiceThread){this_thread::sleep_for (chrono::microseconds(1));}
			pathMutex.lock();
			{
				//~ fwrite((toWrite).c_str(), sizeof(char), toWrite.size(), pathFilef);
				if(compression){
					pathCompressed->write(toWrite.c_str(),toWrite.size());
				}else{
					fwrite((toWrite).c_str(), sizeof(char), toWrite.size(), pathFilef);
				}
				threadToPrint=(threadToPrint+1)%coreNumber;
			}
			pathMutex.unlock();
		}else{
			if(not compressionMode){
				pathMutex.lock();
				{
					if(compression){
						pathCompressed->write(toWrite.c_str(),toWrite.size());
					}else{
						fwrite((toWrite).c_str(), sizeof(char), toWrite.size(), pathFilef);
					}
				}
				pathMutex.unlock();
			}else{
				for(uint i(0);i<nbBuckets;++i){
					pathMutexComp[i].lock();
					{
						fwrite((toWriteComp[i]).c_str(), sizeof(char), toWriteComp[i].size(), pathFileComp[i]);
						toWriteComp[i]="";
					}
					pathMutexComp[i].unlock();
				}
			}
		}
		progressMutex.lock();
		if(++iterLoop%100==0){
			cout<<"Reads: "<<intToString(readNumber)<<endl;
			cout<<"Not anchored : "<<intToString(noOverlapRead)<<" Percent: "<<(100*float(noOverlapRead))/readNumber<<endl;
			cout<<"Anchored and aligned : "<<intToString(alignedRead)<<" Percent: "<<(100*float(alignedRead))/(alignedRead+notAligned)<<endl;
			cout<<"Not aligned: "<<intToString(notAligned)<<" Percent: "<<(100*float(notAligned))/(alignedRead+notAligned)<<endl;
			if(pairedMode){
				cout<<"Super reads: "<<intToString(superReads)<<endl;
				cout<<"Failed super reads compaction: "<<intToString(notCompatedSR)<<endl;
			}
			cout<<endl;
		}else{
		}
		progressMutex.unlock();
	}
}



