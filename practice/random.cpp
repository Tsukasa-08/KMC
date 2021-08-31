#include <random>
#include <iomanip>
#include <iostream>

using namespace std;


int main()
{
	random_device seed_gen;
	mt19937 engine(seed_gen());
	uniform_int_distribution<int> dist(1,64);
	for (int i = 0; i < 5; i++){
		int a = dist(engine);
		int b = dist(engine);

		cout << a << endl;
		cout << b << endl;
	}

/*	random_device seed_gen2;
	mt19937 engine2(seed_gen2());
	uniform_int_distribution<int> dist2(1,64);
	for (int i = 0; i < 5; i++){
		cout << dist2(engine2) << endl;
	}

	random_device rd;
	mt19937 random_engine(rd());
	uniform_real_distribution<double> distr(0,1);

	for (int n = 0; n < 5; ++n) {

		cout << setprecision(6) << distr(random_engine) << endl;
	}
	*/
}

