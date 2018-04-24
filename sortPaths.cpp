#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <math.h>



using namespace std;



struct read_encoding{
	int anchor_number;
	int read_position;
	string path_direction;
	//~ string missmatches_encoding;
};



int charToInt(char c) {
	switch (c) {
		case 'A': return 1;
		case 'C': return 2;
		case 'G': return 3;
	}
	return 4;
}



bool acompare(read_encoding y, read_encoding x) {
	if(abs(x.anchor_number)>abs(y.anchor_number)){return true;}
	if(abs(x.anchor_number)<abs(y.anchor_number)){return false;}

	if(x.read_position>y.read_position){return true;}
	if(x.read_position<y.read_position){return false;}

	if(x.path_direction>y.path_direction){return true;}
	if(x.path_direction<y.path_direction){return false;}



	return false;
}



uint stringToInt(const string& str){
	uint res(0), acc(1);
	for(uint i(0);i<str.size();++i){
		res+=charToInt(str[i])*acc;
		acc*=4;
	}
	return res;
}




int main(int argc, char ** argv){
	if(argc<2){
		cout<<"[path file]"<<endl;
		exit(0);
	}
	string input(argv[1]);

	ifstream in(input);
	if(not in.good()){
		cout<<"Problem with file opening"<<endl;
		exit(1);
	}
	//FILE OPEN
	vector<read_encoding> reads;
	//INPUT

	string parse;
	int anchor_number,read_position;
	while(not in.eof()){
		getline(in,parse,':');
		if(not parse.empty()){
			anchor_number=(stoi(parse));
		}else{
			anchor_number=0;
		}
		getline(in,parse,':');
		string path_direction((parse));
		getline(in,parse,':');
		if(not parse.empty()){
			read_position=(stoi(parse));
		}else{
			read_position=0;
		}
		getline(in,parse,'\n');
		//~ string mismatches_encoding(parse);
		reads.push_back({anchor_number,read_position,path_direction});
	}

	sort(reads.begin(),reads.end(),acompare);
	//INPUT SORTED




	for(uint i(0);i<reads.size();++i){
		cout<<reads[i].anchor_number<<":"<<reads[i].read_position<<":"<<reads[i].path_direction<<endl;
	}
}







