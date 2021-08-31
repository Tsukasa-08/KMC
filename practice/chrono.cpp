#include <iostream>
#include <chrono>

using namespace std;

int main() 
{

	chrono::system_clock::time_point start, end;
	time_t time_stamp;

	start = chrono::system_clock::now();
	
	int sum = 0;
	for (int i = 0; i != 10000; i++) {
		sum += i;	
	}
	cout << "sum = " << sum << endl;

	end = chrono::system_clock::now();

	auto time = end - start;

	time_stamp = chrono::system_clock::to_time_t(start);
	cout << ctime(&time_stamp);

	auto msec = chrono::duration_cast<chrono::nanoseconds>(time).count();
	cout << msec << " nsec" << endl;

	return 0;


}
