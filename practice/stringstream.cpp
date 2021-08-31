#include <iostream>
#include <string>
#include <sstream>
#include <typeinfo>

using namespace std;

int main(){

	string s = "abc def ghi jkl mno";
	stringstream ss;

	cout << "string s type = " << typeid(s).name() << endl;


	ss << s;
	string ret;
	while(!ss.eof()){
		ss >> ret;
		cout << ret << endl;
	}
}

