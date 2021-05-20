#include <iostream>
#include <random>
#include <vector>
#include <algorithm>

using namespace std;

int main(){

	vector<int> v(10, 0.0);
	for (int i = 0; i != v.size(); i++) {
		v[i] = i+1;
		cout << "v[" << i << "] = " << v[i] << endl;
		}

	cout << endl;

	random_device get_rand_dev;
	mt19937 get_rand_mt(get_rand_dev());
	shuffle(v.begin(), v.end(), get_rand_mt );

	for (int i = 0; i != v.size(); i++) {
		cout << "v[" << i << "] = " << v[i] << endl;
		}
	
	cout << endl;

	v.resize(5);

	for (int i = 0; i != v.size(); i++) {
		cout << "v[" << i << "] = " << v[i] << endl;
		}

	cout << endl;

	sort(v.begin(),v.end());

	for (int i = 0; i != v.size(); i++) {
		cout << "v[" << i << "] = " << v[i] << endl;
		}
}
