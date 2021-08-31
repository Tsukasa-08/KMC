#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <typeinfo>

using namespace std;

vector<string> split(const string &str, char sep) 
{
	vector<string> v;
	stringstream ss(str);
	string buffer;
	while(getline(ss, buffer, sep)){
		v.push_back(buffer);
	}

	return v;
}

int main(int argc, char* argv[]) {
	ifstream ifs(argv[1]);

	int line_num = 1;
	int const DIRECT_num = 7; //POSCARの"DIRECT"がある行、これ以降に座標一覧
	string line;
	string x;
	string y;
	string z;
	vector<string> v;


	while(getline(ifs,line)){

		//座標一覧まで(7行目まで)はスキップする
		if (line_num <= DIRECT_num){
			line_num += 1;
			continue;
		}

		//座標一覧以降は座標を抽出し、代入する
		else {
			stringstream ss_line;
			ss_line << line;

			while(!ss_line.eof()){
				string ret;
				ss_line >> ret;
				cout << ret << endl;
			}

			cout << line_num << endl;

			line_num += 1;
		}
	}
}

	
