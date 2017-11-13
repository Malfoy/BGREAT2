#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_set>



using namespace std;



struct read_encoding{
	int anchor_number;
	int read_position;
	int path_direction;
	string missmatches_encoding;
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
	if(x.anchor_number>y.anchor_number){
		return true;
	}
	if(x.anchor_number<y.anchor_number){
		return false;
	}
	//~ if(x.path_direction.size()>y.path_direction.size()){
		//~ return true;
	//~ }
	//~ if(x.path_direction.size()<y.path_direction.size()){
		//~ return false;
	//~ }
	if(x.path_direction>y.path_direction){
		return true;
	}
	if(x.path_direction<y.path_direction){
		return false;
	}
	if(x.read_position>y.read_position){
		return true;
	}
	if(x.read_position<y.read_position){
		return false;
	}
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
	vector<read_encoding> reads;
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
		int path_direction(stringToInt(parse));
		getline(in,parse,':');
		if(not parse.empty()){
			read_position=(stoi(parse));
		}else{
			read_position=0;
		}
		getline(in,parse,(char)255);
		string mismatches_encoding(parse);
		reads.push_back({anchor_number,read_position,path_direction,mismatches_encoding});
	}
	sort(reads.begin(),reads.end(),acompare);
	read_encoding pred(reads[0]);
	cout<<pred.anchor_number<<":"<<pred.path_direction<<":"<<pred.read_position<<":"<<pred.missmatches_encoding<<"\n";
	read_encoding read;
	uint max_diff_anchors(0),max_diff_paths(0),max_diff_pos(0),number_anchors(0),number_path(0),number_pos(0);
	for(uint i(1);i<reads.size();++i){

		bool usePrevious(true);
		read=reads[i];
		//~ if(true){
		if(false){
			cout<<read.anchor_number<<":"<<read.path_direction<<":"<<read.read_position<<":"<<read.missmatches_encoding<<"\n";
		}else{
			int diff_anchors(read.anchor_number-pred.anchor_number);
			if(diff_anchors>max_diff_anchors){
				max_diff_anchors=diff_anchors;
			}
			if(diff_anchors==0){
				//~ cout<<":";
				//~ cout<<(uint8_t)diff_anchors;
			}else{
				cout<<"."<<(uint8_t)diff_anchors;
				number_anchors++;
				usePrevious=false;
			}
			if(usePrevious){
				int diff_paths(read.path_direction-pred.path_direction);
				if(diff_paths>max_diff_paths){
					max_diff_paths=diff_paths;
					//~ cerr<<read.path_direction<<" "<<pred.path_direction<<endl;
				}
				if(diff_paths==0){
					//~ cout<<":";
					//~ cout<<(uint8_t)diff_paths;
				}else{
					cout<<"?"<<(uint32_t)diff_paths;
					number_path++;
					usePrevious=false;
				}
			}else{
				//~ cout<<read.path_direction<<":";
				cout<<(uint32_t)read.path_direction;
				number_path++;
			}
			if(usePrevious){
				//TODO look at all diff, code acording to largest, output bitbector
			int diff_pos(read.read_position-pred.read_position);
			number_pos++;
				if(diff_pos>max_diff_pos){
					max_diff_pos=diff_pos;
				}

				if(diff_pos==0){
					//~ cout<<":";
					if(number_pos%2==0)
						cout<<diff_pos;
				}else{
					if(number_pos%2==0)
						cout<<diff_pos;
				}
			}else{
				//~ cout<<read.read_position<<":";
				cout<<read.read_position;
			}
			//~ cout<<read.missmatches_encoding<<"\n";
			//~ cout<<"\n";
			pred=read;
		}
	}
	cerr<<log2(max_diff_anchors)<<" "<<log2(max_diff_paths)<<" "<<log2(max_diff_pos)<<endl;
	cerr<<number_anchors<<" "<<number_path<<" "<<number_pos<<endl;
}







