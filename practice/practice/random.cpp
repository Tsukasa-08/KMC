#include <random>
#include <iostream>

using namespace std;

int main()
{
	random_device seed_gen;
	mt19937 engine(seed_gen());
	uniform_int_distribution<int> dist(1,64);
	for (int i = 0; i < 5; i++){
		cout << dist(engine) << endl;
	}
	cout << "next" << endl;

	random_device seed_gen2;
	mt19937 engine2(seed_gen2());
	uniform_int_distribution<int> dist2(1,64);
	for (int i = 0; i < 5; i++){
		cout << dist2(engine2) << endl;
	}
}

