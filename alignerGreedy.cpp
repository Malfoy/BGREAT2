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



vector<uNumber> Aligner::alignReadGreedyAnchors(const string& read, uint errorMax,const pair<pair<uint,uint>,uint>& anchor){
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
			uint errors(missmatchNumber(read.substr(positionRead-positionUnitig,unitig.size()),unitig,errorMax));
			//~ cout<<read.substr(positionRead-positionUnitig,unitig.size())<<" "<<unitig<<endl;
			if(errors<=errorMax){
				pathBegin={};
				uint errorBegin;
				errorBegin=(checkBeginGreedy(read,{str2num(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
				if(errorBegin+errors<=errorMax){
					pathEnd={(int)unitigNumber};
					uint errorsEnd;
					errorsEnd=(checkEndGreedy(read,{str2num(unitig.substr(unitig.size()-k+1,k-1)),positionRead-positionUnitig+unitig.size()},pathEnd,errorMax-errors-errorBegin));
					if(errorBegin+errors+errorsEnd<=errorMax){
						reverse(pathBegin.begin(),pathBegin.end());
						pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
						return pathBegin;
					}
				}
			}
		}else{
			//CASE 2 : unitig overap read
			//~ cout<<2<<endl;
			uint errors(missmatchNumber(read.substr(positionRead-positionUnitig),unitig.substr(0,read.size()-positionRead+positionUnitig),errorMax));
			if(errors<=errorMax){
				pathBegin={};
				uint errorBegin;
				errorBegin=(checkBeginGreedy(read,{str2num(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
				if(errorBegin+errors<=errorMax){
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
			uint errors(missmatchNumber(unitig.substr(positionUnitig-positionRead),read.substr(0,unitig.size()+positionRead-positionUnitig),errorMax));
			if(errors<=errorMax){
				pathEnd={(int)positionUnitig-(int)positionRead,(int)unitigNumber};
				uint errorsEnd;
				errorsEnd=(checkEndGreedy(read,{str2num(unitig.substr(unitig.size()-k+1,k-1)),positionRead-positionUnitig+unitig.size()},pathEnd,errorMax-errors));
				if(errors+errorsEnd<=errorMax){
					return pathEnd;
				}
			}
		}else{
			//CASE 4 : read included in unitig
			//~ cout<<"4:"<<endl;
			uint errors(missmatchNumber(unitig.substr(positionUnitig-positionRead,read.size()),read,errorMax));
			if(errors<=errorMax){
				return {(int)positionUnitig-(int)positionRead,(int)unitigNumber};
			}
		}
	}
	return {};
}



vector<uNumber> Aligner::alignReadGreedyAnchorsstr(const string& read, uint errorMax, const pair<pair<uint,uint>,uint>& anchor){
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
			uint errors(missmatchNumber(read.substr(positionRead-positionUnitig,unitig.size()),unitig,errorMax));
			if(errors<=errorMax){
				pathBegin={};
				uint errorBegin;
				errorBegin=(checkBeginGreedy(read,{(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
				if(errorBegin+errors<=errorMax){
					pathEnd={(int)unitigNumber};
					uint errorsEnd;
					errorsEnd=(checkEndGreedy(read,{(unitig.substr(unitig.size()-k+1,k-1)),positionRead-positionUnitig+unitig.size()},pathEnd,errorMax-errors-errorBegin));
					if(errorBegin+errors+errorsEnd<=errorMax){
						//~ ++alignedRead;
						reverse(pathBegin.begin(),pathBegin.end());
						pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
						return pathBegin;
					}
				}
			}
		}else{
			//~ if(positionUnitig+anchorSize<=k-1){
				//~ if(positionUnitig+read.size()-positionRead<=k-1){
					//~ return{};
				//~ }
			//~ }
			uint errors(missmatchNumber(read.substr(positionRead-positionUnitig),unitig.substr(0,read.size()-positionRead+positionUnitig),errorMax));
			if(errors<=errorMax){
				pathBegin={};
				uint errorBegin;
				errorBegin=(checkBeginGreedy(read,{(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
				if(errorBegin+errors<=errorMax){
					//~ ++alignedRead;
					reverse(pathBegin.begin(),pathBegin.end());
					pathBegin.push_back((int)unitigNumber);
					return pathBegin;
				}
			}
		}
	}else{
		if(read.size()-positionRead>=unitig.size()-positionUnitig){
			//~ if(positionUnitig>=unitig.size()-k+1){
				//~ if(positionRead<=positionUnitig){
					//~ return{};
				//~ }
			//~ }
			uint errors(missmatchNumber(unitig.substr(positionUnitig-positionRead),read.substr(0,unitig.size()+positionRead-positionUnitig),errorMax));
			if(errors<=errorMax){
				pathEnd={(int)positionUnitig-(int)positionRead,(int)unitigNumber};
				uint errorsEnd;
				errorsEnd=(checkEndGreedy(read,{unitig.substr(unitig.size()-k+1,k-1),positionRead-positionUnitig+unitig.size()},pathEnd,errorMax-errors));
				if(errors+errorsEnd<=errorMax){
					//~ ++alignedRead;
					return pathEnd;
				}
			}
		}else{
			//~ cout<<4<<endl;
			uint errors(missmatchNumber(unitig.substr(positionUnitig-positionRead,read.size()),read,errorMax));
			if(errors<=errorMax){
				//~ ++alignedRead;
				return {(int)positionUnitig-(int)positionRead,(int)unitigNumber};
			}
		}
	}
	return {};
}



uint Aligner::mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	if(overlap.second<=trimingBases){path.push_back(0);return 0;}
	string unitig,readLeft(read.substr(0,overlap.second)),nextUnitig;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getEnd(overlap.first));
	uint miniMiss(errors+1),miniMissIndice(9);
	bool ended(false);
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
			}else if(miss<=miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnLeftEndGreedy(read , path, {str2num(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)}, errors);

			}else if(miss==miniMiss and noMultiMapping and miniMiss<errors+1){
				return errors+1;
			}else if(miss<miniMiss){
				kmer overlapNum(str2num(unitig.substr(0,k-1)));
				if(miss<miniMiss){
					ended=false;
					miniMiss=miss;
					miniMissIndice=i;
					nextUnitig=unitig;
					nextOverlap=overlapNum;
				}
			}
		}
	}

	if (miniMiss<=errors){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(ended){
			path.push_back(offset);
			return miniMiss;
		}
		return miniMiss+mapOnLeftEndGreedy(read , path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)}, errors-miniMiss);
	}
	return miniMiss;
}



uint Aligner::mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<string, uint>& overlap , uint errors){
	if(overlap.second<=trimingBases){path.push_back(0);return 0;}
	string unitig,readLeft(read.substr(0,overlap.second)),nextUnitig;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getEnd(overlap.first));
	uint miniMiss(errors+1),miniMissIndice(9);
	bool ended(false);
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
			}else if(miss<=miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnLeftEndGreedy(read , path, {(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)}, errors);
			}else if(miss==miniMiss and noMultiMapping and miniMiss<errors+1){
				return errors+1;
			}else if(miss<miniMiss){
				ended=false;
				miniMiss=miss;
				miniMissIndice=i;
				nextUnitig=unitig;
				nextOverlap=((unitig.substr(0,k-1)));
			}
		}
	}
	if (miniMiss<=errors){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(ended){
			path.push_back(offset);
			return miniMiss;
		}
		return miniMiss+mapOnLeftEndGreedy(read , path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)}, errors-miniMiss);
	}
	return miniMiss;
}



uint Aligner::mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	string unitig,readLeft(read.substr(overlap.second)),nextUnitig;
	if(readLeft.size()<trimingBases){return 0;}
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=getBegin(overlap.first);
	uint miniMiss(errors+1), miniMissIndice(9);
	bool ended(false);
	kmer nextOverlap(0);
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return 0;
			}else if(miss<=miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return (mapOnRightEndGreedy(read , path, {str2num(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)}, errors));
			}else if(miss==miniMiss and noMultiMapping and miniMiss<errors+1){
				return errors+1;
			}else if(miss<miniMiss){
				ended=false;
				kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
				miniMiss=miss;
				miniMissIndice=i;
				nextUnitig=unitig;
				nextOverlap=overlapNum;
			}
		}
	}
	if(miniMiss<=errors){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(ended){return miniMiss;}
		return (miniMiss+mapOnRightEndGreedy(read , path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)}, errors-miniMiss));
	}
	return miniMiss;
}



uint Aligner::mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<string, uint>& overlap , uint errors){
	string unitig,readLeft(read.substr(overlap.second)),nextUnitig;
	if(readLeft.size()<trimingBases){return 0;}
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=getBegin(overlap.first);
	uint miniMiss(errors+1), miniMissIndice(9);
	bool ended(false);
	string nextOverlap;
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return 0;
			}else if(miss<=miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return (mapOnRightEndGreedy(read , path, {(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)}, errors));
			}else if(miss==miniMiss and noMultiMapping and miniMiss<errors+1){
				return errors+1;
			}else if(miss<miniMiss){
				ended=false;
				miniMiss=miss;
				miniMissIndice=i;
				nextUnitig=unitig;
				nextOverlap=(unitig.substr(unitig.size()-k+1,k-1));
			}
		}
	}
	if(miniMiss<=errors){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(ended){return miniMiss;}
		return (miniMiss+mapOnRightEndGreedy(read , path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)}, errors-miniMiss));
	}
	return miniMiss;
}



uint Aligner::checkBeginGreedy(const string& read,const pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors){
	if(overlap.second<=trimingBases){path.push_back(0);return 0;}
	string readLeft(read.substr(0,overlap.second)),unitig,nextUnitig;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getEnd(overlap.first));
	uint minMiss(errors+1),indiceMinMiss(9);
	bool ended(false);
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
			}else if(miss<=minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnLeftEndGreedy(read, path, {str2num(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)},errors);
			}else if(miss==minMiss and noMultiMapping and minMiss<errors+1){
				return errors+1;
			}else if(miss<minMiss){
				kmer overlapNum(str2num(unitig.substr(0,k-1)));
				ended=false;
				minMiss=miss;
				indiceMinMiss=i;
				nextUnitig=unitig;
				nextOverlap=overlapNum;
			}

		}
	}
	if(minMiss<=errors){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){
			path.push_back(offset);
			return minMiss;
		}
		return minMiss+mapOnLeftEndGreedy(read, path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)},errors-minMiss);
	}
	return minMiss;
}



uint Aligner::checkBeginGreedy(const string& read,const pair<string, uint>& overlap, vector<uNumber>& path, uint errors){
	if(overlap.second<=trimingBases){path.push_back(0);return 0;}
	string readLeft(read.substr(0,overlap.second)),unitig,nextUnitig;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getEnd(overlap.first));
	uint minMiss(errors+1),indiceMinMiss(9);
	bool ended(false);
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
			}else if(miss<=minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnLeftEndGreedy(read, path, {(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)},errors);
			}else if(miss==minMiss and noMultiMapping and minMiss<errors+1){
				return errors+1;
			}else if(miss<minMiss){
				ended=false;
				minMiss=miss;
				indiceMinMiss=i;
				nextUnitig=unitig;
				nextOverlap=((unitig.substr(0,k-1)));
			}

		}
	}
	if(minMiss<=errors){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){
			path.push_back(offset);
			return minMiss;
		}
		return minMiss+mapOnLeftEndGreedy(read, path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)},errors-minMiss);
	}
	return minMiss;
}



uint Aligner::checkEndGreedy(const string& read, const pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors){
	string readLeft(read.substr(overlap.second)),unitig,nextUnitig;
	if(readLeft.size()<=trimingBases){return 0;}
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getBegin(overlap.first));
	uint minMiss(errors+1),indiceMinMiss(9);
	bool ended(false);
	kmer nextOverlap(0);
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=rangeUnitigs[i].first;
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return 0;
			}else if(miss<=minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnRightEndGreedy(read, path, {str2num(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)},errors);
			}else if(miss==minMiss and noMultiMapping and minMiss<errors+1){
				return errors+1;
			}else if(miss<minMiss){
				kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
				minMiss=miss;
				indiceMinMiss=i;
				nextOverlap=overlapNum;
				nextUnitig=unitig;
				ended=false;
			}
		}
	}
	if(minMiss<=errors){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){return minMiss;}
		return minMiss+mapOnRightEndGreedy(read, path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)},errors-minMiss);
	}
	return minMiss;
}



uint Aligner::checkEndGreedy(const string& read, const pair<string, uint>& overlap, vector<uNumber>& path, uint errors){
	string readLeft(read.substr(overlap.second)),unitig,nextUnitig;
	if(readLeft.size()<=trimingBases){return 0;}
	vector<pair<string,uNumber>> rangeUnitigs(getBegin(overlap.first));
	uint minMiss(errors+1),indiceMinMiss(9);
	bool ended(false);
	string nextOverlap;
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=rangeUnitigs[i].first;
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return 0;
			}else if(miss<=minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnRightEndGreedy(read, path, {(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)},errors);
			}else if(miss==minMiss and noMultiMapping and minMiss<errors+1){
				return errors+1;
			}else if(miss<minMiss){
				nextOverlap=((unitig.substr(unitig.size()-k+1,k-1)));
				minMiss=miss;
				indiceMinMiss=i;
				nextUnitig=unitig;
				ended=false;
			}
		}
	}
	if(minMiss<=errors){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){return minMiss;}
		return minMiss+mapOnRightEndGreedy(read, path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)},errors-minMiss);
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
	if(listAnchors.empty()){
		++noOverlapRead;
		return;
	}
	//~ random_shuffle ( listAnchors.begin(), listAnchors.end() );
	while(errors<=errorsMax){
		bool found(false);
		pathMem={};
		for(uint i(0);i<listAnchors.size();++i){
			path={};
			if(stringMode){
				path=alignReadGreedyAnchorsstr(read,errors,listAnchors[i]);
			}else{
				path=alignReadGreedyAnchors(read,errors,listAnchors[i]);
			}
			//MAPPING IS FOUND
			if(not path.empty()){
				if(noMultiMapping){
					if(found){
						pathMem=inclued(path,pathMem);
						if(pathMem.empty()){path={};return;}
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
		}
		++errors;
		if(perfect or found){
			break;
		}
	}
	path=pathMem;
	if(not path.empty()){
		++alignedRead;
	}
}



void Aligner::alignReadFrom(const string& read, vector<int>& path, int unumber){
	vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	for(uint i(0);i<listAnchors.size();++i){
		auto  anchor=listAnchors[i];
		if(anchor.first.first==unumber or anchor.first.first==-unumber){
			path={};
			if(stringMode){
				path=alignReadGreedyAnchorsstr(read,0,anchor);
			}else{
				path=alignReadGreedyAnchors(read,0,anchor);
			}
			return;
		}
	}
}




void Aligner::alignPartGreedy(uint indiceThread){
	vector<pair<string,string>> multiread;
	vector<uNumber> path,path2;
	string read,read2,header,header2,corrected,superRead,toWrite;
	pair<string,string> superpath;
	while(!readFile.eof()){
		toWrite="";
		if(keepOrder){
			while(threadToRead!=indiceThread){
				this_thread::sleep_for (chrono::microseconds(1));
			}
			readMutex.lock();
			getReads(multiread,10000);
			threadToRead=(threadToRead+1)%coreNumber;
			readMutex.unlock();
		}else{
			readMutex.lock();
			getReads(multiread,10000);
			readMutex.unlock();
		}
		if(pairedMode){
			for(uint i(0);i+1<multiread.size();i+=2){
				header=multiread[i].first;
				read=multiread[i].second;
				header2=multiread[i+1].first;
				read2=multiread[i+1].second;
				readNumber+=2;
				bool rc(false), noOverlap(false), overlapFound(false);
				alignReadOpti(read,path);
				alignReadOpti(read2,path2);
				if(path.empty()){++notAligned;}
				if(path2.empty()){++notAligned;}
				superpath=(recoverSuperReadsPairedNoStr(path,path2));
				if(superpath.first!=""){
					if(superpath.second!=""){
						toWrite+=header+'\n'+superpath.first+'\n'+header2+'\n'+superpath.second+'\n';
					}else{
						toWrite+=header+'\n'+superpath.first+'\n';
					}
				}else{
					if(superpath.second!=""){
						toWrite+=header2+'\n'+superpath.second+'\n';
					}
				}
			}
		}else{
			for(uint i(0);i<multiread.size();i++){
				//~ cout<<"go"<<endl;
				header=multiread[i].first;
				read=multiread[i].second;
				++readNumber;
				alignReadOpti(read,path);
				if(path.empty()){++notAligned;}
				if(correctionMode){
					if(not path.empty()){
						uint position(path[0]);
						path=vector<uNumber>(&path[1],&path[path.size()]);
						superRead=(recoverSuperReads(path));
						if(superRead.substr(position,read.size())!=read){
							cout<<superRead.substr(position,read.size())<<endl;
							cout<<read<<endl;
							cout<<"NOOOO"<<endl;
							cin.get();
						}
						if(superRead!=""){
							toWrite+=header+'\n'+superRead.substr(position,read.size())+'\n';
						}else{
							toWrite+=header+'\n'+read+'\n';
						}
					}else{
						toWrite+=header+'\n'+read+'\n';
					}
				}else{
					if(not path.empty()){
						if(preciseOutput){
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
							superRead=(recoverSuperReadsNoStr(path));
							if(superRead!=""){
								toWrite+=header+'\n'+superRead+'\n';
							}
						}else{
							path=vector<uNumber>(&path[1],&path[path.size()]);
							superRead=(recoverSuperReadsNoStr(path));
							if(superRead!=""){
								toWrite+=header+'\n'+superRead+'\n';
							}
						}
					}else{
						//~ cerr<<header<<endl;
						//~ cerr<<read<<endl;
					}
				}
			}
		}
		if(keepOrder){
			while(threadToPrint!=indiceThread){this_thread::sleep_for (chrono::microseconds(1));}
			pathMutex.lock();
			{
				fwrite((toWrite).c_str(), sizeof(char), toWrite.size(), pathFilef);
				threadToPrint=(threadToPrint+1)%coreNumber;
			}
			pathMutex.unlock();
		}else{
			pathMutex.lock();
			{
				fwrite((toWrite).c_str(), sizeof(char), toWrite.size(), pathFilef);
			}
			pathMutex.unlock();
		}
	}
}



