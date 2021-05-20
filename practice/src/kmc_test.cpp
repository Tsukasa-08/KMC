#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <algorithm>

using namespace std;

//サイトクラスを定義
class Site {

private:
	int site_id;
	int site_atom;
	vector<double> site_frac_coods;

public:
	
	//デフォルトコンストラクタ デフォルトでは空孔にしておく
	Site() : site_id(0), site_frac_coods(3,0.0), site_atom(1){
	}

	//セッタ
	void set_site_id(int id) { site_id = id; };

	void set_site_atom(int atom) { site_atom = atom; };

	void set_site_frac_coods(vector<double> v) { site_frac_coods = v; };

	//ゲッタ
	int get_site_id() { return site_id; } ;

	int get_site_atom() { return site_atom; } ;

	vector<double> get_site_frac_coods() { return site_frac_coods; } ;
};


//ジャンプクラスを定義
class Jump {

private:
	int start_site_id;
	int start_site_atom;
	int end_site_id;
	int end_site_atom;
	double freq;

public:
	//コンストラクタ
	Jump(int s_s_id = 0, int s_s_atom = 0, int e_s_id = 0, int e_s_atom = 0, double fq = 0){

	};


	//セッタ	
	void set_start_site_id(int id) { start_site_id = id; };

	void set_start_site_atom(int atom){ start_site_atom = atom; };

	void set_end_site_id(int id){ end_site_id = id; };

	void set_end_site_atom(int atom){ end_site_atom = atom; };

	void set_freq(double fq){ freq = fq; };

	//ゲッタ
	int get_start_site_id() { return start_site_id; }

	int get_start_site_atom() { return start_site_atom; }

	int get_end_site_id() { return end_site_id; }

	int get_end_site_atom() { return end_site_atom; }

	double get_freq() { return freq; }
};



//プロトンクラスを定義
class Proton {

private:
	int proton_id;
	int proton_start_site;
	int proton_end_site;
	int cnt_periodicboundary;
};




//ファイルを読み込むためのsplit関数を定義

vector<string> split(const string& input, char delimiter)
{
	istringstream stream(input);
	string field;
	vector<string> result;
	while (getline(stream, field, delimiter)){
		result.push_back(field);
	}
	return result;

}


//原子種の数字を定義

int const number_vacancy = 1;
int const number_proton = 2;

//メイン関数の開始
int main()
{

	//JMPDATAを読み込む
	ifstream ifs1("JMPDATA");

	int	const jump_total_number = 256;
	vector<Jump> jump_vectors(jump_total_number, Jump());
	string line;
	int n_lines = 0;
	while (getline(ifs1, line)){


		vector<string> strvec = split(line, ',');

		jump_vectors[n_lines].set_start_site_id(stoi(strvec[0]));
		jump_vectors[n_lines].set_start_site_atom(stoi(strvec[1]));
		jump_vectors[n_lines].set_end_site_id(stoi(strvec[2]));
		jump_vectors[n_lines].set_end_site_atom(stoi(strvec[3]));
		jump_vectors[n_lines].set_freq(stod(strvec[4]));

//		cout << jump_vectors[n_lines].get_freq() << endl;
		
		n_lines += 1;
	}	
			

	//POSCARを読み込む
	ifstream ifs2("POSCAR");

	n_lines = 1;
	int const DIRECT_num = 7; //POSCARの"DIRECT"がある行、これ以降に座標一覧がある
	int const site_total_number = 64;
	vector<Site> sites(site_total_number,Site());
		
	//POSCARを1行ごとに読み込んでいく
	while(getline(ifs2,line)){

		//座標一覧までは飛ばす
		if(n_lines <= DIRECT_num){
			n_lines += 1;
			continue;
		}

		//座標一覧に到達以降
		else {
			stringstream ss_line; 
			ss_line << line;

			//1行を空白で3つの座標に分割していく
			while(!ss_line.eof()){

				//座標をコピーするためのfrac_dblvecを作成する
				vector<string> frac_strvec(3,"example");
				vector<double> frac_dblvec(3,0.0);
				
				for (int i=0; i <= 2; i++) {
					
					ss_line >> frac_strvec[i];
					frac_dblvec[i] = stod(frac_strvec[i]);

				//	cout << frac_dblvec[i] << endl;
				}	

				//siteのidは「現在の行数-"DIRECT"の行数」
				int site_id_tmp = n_lines - DIRECT_num;

				//sites自体は0から始まるので、1つずらして代入する
				sites[site_id_tmp-1].set_site_id(site_id_tmp);
				sites[site_id_tmp-1].set_site_frac_coods(frac_dblvec);

			//	cout << sites[site_id_tmp-1].get_site_id() << endl;
			//	for (int i=0; i <3; i++) {
			//		cout << sites[site_id_tmp-1].get_site_frac_coods()[i] << endl;
				}
			
		}

		n_lines += 1;
		
	}

	//初期配置を生成する
	//まずはプロトンを配置する数を決定する
	random_device rnd;
	mt19937 mt(rnd());
	uniform_int_distribution<int> rnd_p_place_n(1,site_total_number);

	int p_place_n = rnd_p_place_n(mt);
	cout << "プロトン配置数" << p_place_n << endl;

	//次に、どのサイトにプロトンを配置するかを決める
	//サイトの数を要素にもつvectorを生成
	vector<int> rnd_p_place_v(site_total_number, 0.0);
	for (int i = 0; i != rnd_p_place_v.size(); i++) {
		rnd_p_place_v[i] = i + 1;
	}

	//乱数を生成し、vectorをシャッフルする
	random_device get_rand_dev;
	mt19937 get_rand_mt(get_rand_dev());
	shuffle( rnd_p_place_v.begin(), rnd_p_place_v.end(), get_rand_mt );
	
	//先頭からプロトンを配置する数分(p_place_n)だけ抜き出しソートする
	rnd_p_place_v.resize(p_place_n);
	sort( rnd_p_place_v.begin(), rnd_p_place_v.end() );

	for (int i = 0; i != rnd_p_place_v.size(); i++) {
		cout << "rnd_p_place_v[" << i << "] = " << rnd_p_place_v[i] << endl;
	}

}


