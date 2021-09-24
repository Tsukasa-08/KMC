#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <regex>

#include "Site.h"
#include "Jump.h"
#include "Diffusionspecie.h"
#include "param.hpp"

using namespace std;

//
/*
//サイトクラスを定義
class Site {

private:
	int site_id;
	int site_atom;
	int diffusion_id;
	vector<double> site_frac_coords;

public:
	
	//デフォルトコンストラクタ デフォルトでは空孔にしておく
	Site() : site_id(0), site_frac_coords(3,0.0), site_atom(1), diffusion_id(-1){
	}

	//セッタ
	void set_site_id(int id) { site_id = id; };

	void set_site_atom(int atom) { site_atom = atom; };

	void set_diffusion_id(int d_id) { diffusion_id = d_id; };

	void set_site_frac_coords(vector<double> v) { site_frac_coords = v; };

	//ゲッタ
	int get_site_id() { return site_id; } ;

	int get_site_atom() { return site_atom; } ;

	int get_diffusion_id() { return diffusion_id; } ;

	vector<double> get_site_frac_coords() { return site_frac_coords; } ;
};


//ジャンプクラスを定義
class Jump {

private:
	int start_site_id;
	int start_site_atom;
	int end_site_id;
	int end_site_atom;
	vector<double> jump_vector;
	double freq;

public:
	//デフォルトコンストラクタ
	Jump() : start_site_id(0), start_site_atom(0), end_site_id(0), end_site_atom(0),jump_vector(3,0.0), freq(0.0){

	};


	//セッタ	
	void set_start_site_id(int id) { start_site_id = id; }

	void set_start_site_atom(int atom){ start_site_atom = atom; }

	void set_end_site_id(int id){ end_site_id = id; }

	void set_end_site_atom(int atom){ end_site_atom = atom; }

	void set_jump_vector(vector<double> jv) { jump_vector = jv; }

	void set_freq(double fq){ freq = fq; }

	//ゲッタ
	int get_start_site_id() { return start_site_id; }

	int get_start_site_atom() { return start_site_atom; }

	int get_end_site_id() { return end_site_id; }

	int get_end_site_atom() { return end_site_atom; }

	vector<double> get_jump_vector() { return jump_vector; }

	double get_freq() { return freq; }
};



//拡散種クラスを定義
class Diffusionspecie {

private:
	int diffusion_id;
	vector<double> jump_total;

public:

	Diffusionspecie() : diffusion_id(-1), jump_total(3,0.0) {
	};

	//セッタ
	void set_diffusion_id(int d_id) { diffusion_id = d_id; }

	void set_jump_total(vector<double> jt) { jump_total = jt; }

	//ゲッタ
	int get_diffusion_id() { return diffusion_id; }

	vector<double> get_jump_total() { return jump_total; }
};

*/





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

//分率座標をcartesian座標に変換するtranscoords関数を定義

Eigen::Vector3d transcoords(const vector<double> vector_frac, Eigen::Matrix3d lattice_matrix)
{
	Eigen::Vector3d vector_frac_eigen;
	vector_frac_eigen << vector_frac[0], vector_frac[1], vector_frac[2];
	Eigen::Vector3d vector_cartesian;
	vector_cartesian = lattice_matrix*vector_frac_eigen;

	return vector_cartesian;
}



//原子種の数字を定義

int const number_vacancy = 1;
int const number_proton = 2;


//物理定数を定義

double const kb = 1.380649 * pow(10,-23);
double const ec = 1.60217662 * pow(10,-19);

//プロトンの電荷を定義、要改善
double q_charge = 1 * ec;


//メイン関数の開始
int main()
{

	//実行時間を計測する
	chrono::system_clock::time_point start, end;
	time_t time_stamp;


	start = chrono::system_clock::now();

	//make DiffusionCoefficient vector
	vector< vector<double> > D_t_3d_vector;
	vector< vector<double> > D_j_3d_vector;
	vector< vector<double> > D_c_3d_vector;

	//make ElectricalConductivity vector
	vector< vector<double> > Sigma_x_vector;

	//INPUTを読み込む
	param::parameter param("INPUT");
	int step_max = param.get<int>("NSTEPS", 0);
	int loop_max = param.get<int>("NLOOPS", 0);
	int p_place_n = param.get<int>("NDIFFS", 0);
	int E_field_yes = param.get<int>("EFIELDON", 0);
	double E_field_strength_for_pow = param.get<double>("EFIELD", 0);
	double E_field_strength = pow(10, E_field_strength_for_pow);
	
	double temperture = param.get<double>("TEMP", 0);
	int E_field_axis = param.get<int>("AXIS", 0);
	
	//入力確認用logファイルを作成
	ofstream ofs_log("log", ios::app);

	//読み込めたか確認用
	if (!param) {
		cerr << "Could not find file INPUT" << endl;
		abort();
	}
	cout << "INPUTファイルを読み込み" << endl;
	cout << "\t" << "KMCのステップ数(何回行うか) NSTEPS =" << step_max << endl;
	cout << "\t" << "KMCのループ数(1回のKMCで何回イベントを起こすか) NLOOPS =" << loop_max << endl;
	cout << "\t" << "拡散種をいくつ配置するか NDIFFS = " << p_place_n << endl;
	cout << "\t" << "電場の強さ EFIELD = " << E_field_strength << endl;
	cout << "\t" << "温度 TEMP = " << temperture << endl;
	cout << "\t" << "電場方向 " << E_field_axis << " (x=0, y=1, z=2)" << endl;
	cout << endl;

	ofs_log << "INPUTファイルを読み込み" << endl;
	ofs_log << "\t" << "KMCのステップ数(何回行うか) NSTEPS =" << step_max << endl;
	ofs_log << "\t" << "KMCのループ数(1回のKMCで何回イベントを起こすか) NLOOPS =" << loop_max << endl;
	ofs_log << "\t" << "拡散種をいくつ配置するか NDIFFS = " << p_place_n << endl;
	ofs_log << "\t" << "電場の強さ EFIELD = " << E_field_strength << endl;
	ofs_log << "\t" << "温度 TEMP = " << temperture << endl;
	ofs_log << "\t" << "電場方向 " << E_field_axis << " (x=0, y=1, z=2)" << endl;
	ofs_log << endl;

	
	//JMPDATAを読み込む
	//まずは行数を取得
	ifstream for_line1("JMPDATA");

	//読み込めたか確認用
	if (!for_line1) {
		cerr << "Could not find file JMPDATA" << endl;
		abort();
	}

	int jump_total_number = 0;
	string for_line_reader1;
	while (getline(for_line1,for_line_reader1)) {

		jump_total_number++;
	}

	cout << "JMPDATA read" << endl;
	cout << "\t" << "jump_total_number = " << jump_total_number << endl;
	cout << endl;

	ofs_log << "JMPDATA read" << endl;
	ofs_log << "\t" << "jump_total_number = " << jump_total_number << endl;
	ofs_log << endl;

	//int	const jump_total_number = 256;
	ifstream ifs1("JMPDATA");
	vector<Jump> jumps(jump_total_number, Jump());
	string line;
	int n_lines = 0;
	while (getline(ifs1, line)){


		vector<string> strvec = split(line, ',');

		jumps[n_lines].set_start_site_id(stoi(strvec[0]));
		jumps[n_lines].set_start_site_atom(stoi(strvec[1]));
		jumps[n_lines].set_end_site_id(stoi(strvec[2]));
		jumps[n_lines].set_end_site_atom(stoi(strvec[3]));
		jumps[n_lines].set_freq(stod(strvec[4]));

	//	cout << "jumps[" << n_lines << "].freq = " << jumps[n_lines].get_freq() << endl;
		
		n_lines += 1;
	}	
			

	//POSCARを読み込む

	n_lines = 1;
	
	//まずは行数の情報を取得
	ifstream for_line2("POSCAR");
	
	//読み込めたか確認用
	if (!for_line2) {
		cerr << "Could not find file POSCAR" << endl;
		abort();
	}

	int site_total_number = 0;
	int DIRECT_num = 0;
	string for_line_reader2;
	regex re("direct", regex_constants::icase);
	smatch m;
	while (getline(for_line2,for_line_reader2)) {

		site_total_number++;
		
		//Directがある行を特定しておく、だいたい7行目
		if (regex_search(for_line_reader2, m, re)) {
			DIRECT_num = site_total_number;
			cout << "Directは" << site_total_number << "行目にあります" << endl;
		}
	}

	//サイト数=POSCARの行数-(座標一覧までの行数)
	site_total_number -= DIRECT_num;
	cout << "POSCAR read" << endl;
	cout << "\t" << "site_total_number = " << site_total_number << endl;
	cout << endl;

	ofs_log << "POSCAR read" << endl;
	ofs_log << "\t" << "site_total_number = " << site_total_number << endl;
	ofs_log << endl;

	//POSCARを1行ごとに読み込んでいく
	vector<Site> sites(site_total_number,Site());
	Eigen::Matrix3d lattice_matrix;
	int n_lattice_line = 0;
	
	ifstream ifs2("POSCAR");
	while(getline(ifs2,line)){


		//座標一覧までは飛ばす
		if(n_lines <= DIRECT_num){

			//格子ベクトルのある3~5行目だった場合は読み込む
			if(n_lines >= 3 && n_lines <= 5) {
				
			//	cout << "n_lines = " << n_lines << endl;
				stringstream ss_line;
				ss_line << line;

					for (int i=0; i <= 2; i++) {
						
						string tmp;
						ss_line >> tmp;
						lattice_matrix(n_lattice_line,i) = stod(tmp);

					}	
				
				n_lattice_line += 1;

			}

			n_lines += 1;
			continue;
		}



		//座標一覧に到達以降
		else {
		//	cout << "n_lines = "  << n_lines << endl;

			stringstream ss_line; 
			ss_line << line;
			string s;
			int s_counter = 0;

			//座標をコピーするためのfrac_dblvecを作成する
			vector<double> frac_dblvec(3,0.0);
			//cout << "s = " << s << endl;
			//cout << "s_counter = " << s_counter << endl;

			//1行を空白で3つの座標に分割していく
			while(getline(ss_line, s, ' ')){

				frac_dblvec[s_counter] = stod(s);
			//	cout << "frac_dblvec[" << s_counter << "] = "<< frac_dblvec[s_counter] << endl;
				s_counter += 1;
			}

			//siteのidは「現在の行数-"DIRECT"の行数」
			int site_id_tmp = n_lines - DIRECT_num;

			//sites自体は0から始まるので、1つずらして代入する
			sites[site_id_tmp-1].set_site_id(site_id_tmp);
			sites[site_id_tmp-1].set_site_frac_coords(frac_dblvec);

				//cout << sites[site_id_tmp-1].get_site_id() << endl;
				//for (int i=0; i <3; i++) {
					//cout << sites[site_id_tmp-1].get_site_frac_coords()[i] << endl;
				//}
			
			
		}

		n_lines += 1;
		
	}
	
		//格子ベクトル確認用
		cout << "lattice_matrix = " << endl;
		for (int i = 0; i <= 2; i++) {
			cout << "\t"  << lattice_matrix.row(i) << endl;
		} 
		cout << endl;

		ofs_log << "lattice_matrix = " << endl;
		for (int i = 0; i <= 2; i++) {
			ofs_log << "\t"  << lattice_matrix.row(i) << endl;
		} 
		ofs_log << endl;

/*	//Site確認用
	for (int i = 0; i != sites.size() ; i++) {
		cout << "sites[" << i << "].id = " << sites[i].get_site_id() << endl;
		for (int j = 0; j <= 2; j++){
			cout << '\t' << sites[i].get_site_frac_coords()[j] << endl;
		}
	}
*/





	//Site.site_frac_coordsをもとに、Jump.jump_vectorを生成する
	for (int i = 0; i != jumps.size(); i++) {

		cout << "jumps_count = " << i << endl;
		//始点と終点のsite_idを取得
		int start_site_id_tmp = jumps[i].get_start_site_id();
		//cout << "start_site_id_tmp = " << start_site_id_tmp << endl;
		int end_site_id_tmp = jumps[i].get_end_site_id();
		//cout << "end_site_id_tmp = " << end_site_id_tmp << endl;

		//始点と終点の分率座標を取得
		//sites[i]のsite_idはi+1なので、tmpから1を引いておく
		vector<double> start_frac_coords = sites[start_site_id_tmp-1].get_site_frac_coords();
		//cout << "start_frac_coords = " << start_frac_coords[0] << endl;
		vector<double> end_frac_coords = sites[end_site_id_tmp-1].get_site_frac_coords();
		//cout << "end_frac_coords = " << end_frac_coords[0] << endl;

		//仮のjump_vector_tmpを作成する
		vector<double> jump_vector_tmp(3,0.0);

		//成分ごとに計算していく
		for (int k = 0; k != jump_vector_tmp.size(); k++) {
			jump_vector_tmp[k] = end_frac_coords[k] - start_frac_coords[k];

			//周期的境界条件
			//ベクトルの成分の絶対値が0.5を超えている(=ユニットセルの半分を移動している)場合に
			//1を足し引きして補正する
			if (jump_vector_tmp[k] > 0.5) {
				jump_vector_tmp[k] -= 1;
			}
			else {
				if (jump_vector_tmp[k] < -0.5) {
					jump_vector_tmp[k] += 1;
				}
				else {
					//continue;
				}
			}


		//cout << "jump_vector_tmp[" << k << "] = " << jump_vector_tmp[k] << endl;
		

		//jump_vector_tmpをjump_vectorに代入する
		jumps[i].set_jump_vector(jump_vector_tmp);

/*		//確認用
		for (int j = 0; j != jumps[i].get_jump_vector().size(); j++) {
			cout << "jump[" << i << "].jump_vector[" << j << "] = " << jumps[i].get_jump_vector()[j] << endl;
		}
		*/
		

		}


	}



	//電場がかかっていた場合、ジャンプ頻度を補正する
	if (E_field_yes) {

		//cartesian座標軸方向の単位ベクトルを作成し、電場の大きさをかけて電場ベクトルとする
		//変数E_field_axisによって作成する単位ベクトルを変える(x=0, y=1, z=2)
		
		Eigen::Vector3d E_field_vector(3); 

		switch (E_field_axis) {
			case 0 : E_field_vector << 1,0,0 ; break;
			case 1 : E_field_vector << 0,1,0 ; break;
			case 2 : E_field_vector << 0,0,1 ; break;
		}

		E_field_vector *= E_field_strength;
		cout << "E_field_vector = " << E_field_vector << endl;
		ofs_log << "E_field_vector = " << endl;
		ofs_log <<  E_field_vector << endl;

		//生成したジャンプに対し操作を行っていく
		for (int i = 0; i != jumps.size(); i++) {
			
			//ジャンプベクトルをfracからcartesianに直す(関数作ったので必要なし)
			/*
			Eigen::Vector3d jump_vector_frac;
			Eigen::Vector3d jump_vector_cartesian;
			jump_vector_frac << jumps[i].get_jump_vector()[0], jumps[i].get_jump_vector()[1], jumps[i].get_jump_vector()[2];
			jump_vector_cartesian = lattice_matrix*jump_vector_frac;
			cout << "jump_vector_cartesian = " << jump_vector_cartesian << endl;
			
			*/

			//transcoords関数でジャンプベクトルをfracからcartesianに直す
			Eigen::Vector3d jump_vector_cartesian = transcoords(jumps[i].get_jump_vector(), lattice_matrix);
			
			//ジャンプベクトル方向の電場の大きさを計算する
			double E_j_dot = E_field_vector.dot(jump_vector_cartesian);
			//cout << "E_j_dot = " << E_j_dot << endl;
			double E_along_jump_strength = E_j_dot / jump_vector_cartesian.norm();
			//cout << "E_along_jump_strength = " << E_along_jump_strength << endl;
		
			//ΔEmigを求める、1*で良いのはプロトンのみなことに注意
			double delta_E_mig = q_charge * E_along_jump_strength * jump_vector_cartesian.norm()/2;
			//cout << "delta_E_mig [J] = " << delta_E_mig << endl;

			//ジャンプ頻度を補正する
			double fixed_jump_freq = jumps[i].get_freq() * exp(delta_E_mig/(kb*temperture));
			//cout << "jump_frep = " << jumps[i].get_freq() << endl;
			//cout << "fixed_jump_freq = " << fixed_jump_freq << endl;
			jumps[i].set_freq(fixed_jump_freq);
		}

	}

	else {
		cout << '\t' << "no E_field" << endl;	
		ofs_log << '\t' << "no E_field" << endl;	
	}


	//拡散係数記録用のファイルOUTPUTを作成しておく(必要なくなった)
	//ofstream ofs("OUTPUT", ios::app);


	for (int step_counter = 1; step_counter <= step_max; step_counter++) {
		
		//実行時間を計測する
	//	chrono::system_clock::time_point start, end;
	//	time_t time_stamp;


	//	start = chrono::system_clock::now();



		//確認用
//		cout << "step_counter = " << step_counter << " start " << endl;

		//乱数を準備しておく
		random_device rnd;
		mt19937 mt(rnd());

		//初期配置を生成する

		//まずはプロトンを配置する数を決定する(input.cfgで入力済みのはず)
		//もしinput.cfgでプロトンの配置数が設定されておらず、デフォルトの0が採用されていた場合
		if (p_place_n == 0) {

			cout << "NDIFFS not defined. decide NDIFFS." << endl;
			
			//1から(site_total_number-1)のうち、一様ランダムに1つを決定=プロトン配置数
			uniform_int_distribution<int> rnd_p_place_n(1,site_total_number-1);
			p_place_n = rnd_p_place_n(mt);
			cout << "プロトン配置数" << p_place_n << endl;
		

		}

		//cout << "NDIFFS = " << p_place_n << endl;

		//次に、どのサイトにプロトンを配置するかを決める
		//(サイトの数-1)を要素にもつvectorを生成
		vector<int> rnd_p_place_v(site_total_number-1, 0.0);
		for (int i = 0; i != rnd_p_place_v.size(); i++) {
			rnd_p_place_v[i] = i + 1;
		}

		//乱数を生成し、vectorをシャッフルする
		shuffle( rnd_p_place_v.begin(), rnd_p_place_v.end(), mt );
		
		//先頭からプロトンを配置する数分(p_place_n)だけ抜き出しソートする
		rnd_p_place_v.resize(p_place_n);
		sort( rnd_p_place_v.begin(), rnd_p_place_v.end() );

	/*	//確認用
		for (int i = 0; i != rnd_p_place_v.size(); i++) {
			cout << "rnd_p_place_v[" << i << "] = " << rnd_p_place_v[i] << endl;
		}
	*/	


		//Diffusionspecieクラスのvectorをつくる(数はプロトン配置数=p_place_n)
		vector<Diffusionspecie> diffusion_species(p_place_n);
		//それぞれのidを[i]にたいしてi+1で設定する(1からp_place_nまでdiffusion_idとして通し番号をふる)
		for (int i = 0; i != diffusion_species.size(); i++) {
			diffusion_species[i].set_diffusion_id(i+1);
		}
		
		


		//プロトンを配置するSiteクラスのsite_atomをnumber_proton(今回は2を割当)に変更する
		//と同時に、diffusion_idを設定する
		
		int d_id = 1;

		for (int i = 0; i != rnd_p_place_v.size(); i++) {
			for (int k = 0; k != sites.size(); k++) {

				//sites[k]=Siteクラスのsite_idが、プロトンを配置するsite_idかどうかを判定
				if (sites[k].get_site_id() == rnd_p_place_v[i]) {
					sites[k].set_site_atom(number_proton);
					sites[k].set_diffusion_id(d_id);
					d_id++;

					//cout << "proton at site id = " << sites[k].get_site_id() << endl;
					break;
				}
				else {
					continue;
				}
			}
		}

	//	if ( d_id != p_place_n+1 )
			//cout << "d_id != p_place_n+1" << endl;
	//	else
			//cout << "d_id = p_place_n+1" << endl;

/*		//確認用
		for (int i = 0; i != sites.size(); i++) {
			cout << "site[" << i << "] ";
			cout << "site_id = " << sites[i].get_site_id() << " " ;
			cout << "site_atom = " << sites[i].get_site_atom() << " " ;
			cout << "diffusion_id = " << sites[i].get_diffusion_id() << endl;
		}
*/
		
		
		//シミュレーション時間total_timeを定義
		double total_time = 0.0;

		cout << endl;
		//cout << "\t" << "total_time = " << total_time << endl;


		//メインループの開始
		cout << endl;
		//cout << "\t" << "main loop start" << endl;
		for (int loop_counter = 1; loop_counter <= loop_max; loop_counter++) { 			


			//実行時間を計測する
	//		chrono::system_clock::time_point start, end;
	//		time_t time_stamp;

	//		start = chrono::system_clock::now();


			//確認用
			//cout << endl;
			//cout << "\t" << "loop_counter = " << loop_counter << endl;

			//系で起きうる事象(今回はジャンプ)を列挙し、jumps_possibleに入れていく
			//まずは始点のidと原子種が一致するものを抽出しjumps_possible_tmpに入れる
			vector<Jump> jumps_possible_tmp;
			vector<Jump> jumps_possible;
			for (int i = 0; i != jumps.size(); i++) {
				for (int k = 0; k != sites.size(); k++) {

					if ( sites[k].get_site_id() == jumps[i].get_start_site_id() 
						 && 
						 sites[k].get_site_atom() == jumps[i].get_start_site_atom() )
					{

						jumps_possible_tmp.push_back(jumps[i]);
				//		 cout << "sites[" << k << "]_id = " << sites[k].get_site_id() << "  " ;
				//		 cout << "jumps[" << i << "]_id = " << jumps[i].get_start_site_id() << endl;

					}

				}
			}

			//次に終点のidと原子種が一致するものを抽出しjumps_possibleに入れる
			for (int i = 0; i != jumps_possible_tmp.size(); i++) {
				for (int k = 0; k != sites.size(); k++) {

					if ( sites[k].get_site_id() == jumps_possible_tmp[i].get_end_site_id() 
						 && 
						 sites[k].get_site_atom() == jumps_possible_tmp[i].get_end_site_atom() ) {

						jumps_possible.push_back(jumps_possible_tmp[i]);

					}

				}
			}
		

			//起きうるジャンプについて、freqの総和をとりfreq_sumに格納する
			
			double freq_sum = 0.0;
			for (int i = 0; i != jumps_possible.size(); i++) {
				freq_sum += jumps_possible[i].get_freq();
			}

			//確認用
	/*		for (int i = 0; i != jumps_possible.size(); i++) {
				cout << "\t" << "\t" << "jumps_possible[" << i << "] = " << jumps_possible[i].get_freq() << endl;
			}

			cout << "\t"  << "\t" << "freq_sum = " <<  freq_sum << endl;
			
			*/
			







			//実際に起こるイベントを決める
			//0から1の乱数を生成する
			uniform_real_distribution<double> random0to1(0,1);
			double rho_1 = random0to1(mt);
			//cout << "\t" << "\t" << "rho_1 = " << rho_1 << endl;
			
			//jumps_possibleのfreqを順に足し上げていき、和がfreq_sumを超えたときの整数lを取得する
			double freq_sum_tmp = 0.0;
			int jump_happen_number;
			for (int l = 0; l != jumps_possible.size(); l++) {
				freq_sum_tmp += jumps_possible[l].get_freq();

				if (freq_sum_tmp > rho_1 * freq_sum) {
					jump_happen_number = l;
					break;
				}
			}

			//cout << "\t" << "\t"<< "jump_happen_number = " << jump_happen_number  << endl;



			//実際にイベントを起こす
			int jumps_start_id = jumps_possible[jump_happen_number].get_start_site_id();
			int jumps_end_id = jumps_possible[jump_happen_number].get_end_site_id();
			
			//startとendのsite_atomを交換
			int atom_tmp = sites[jumps_start_id-1].get_site_atom();
			sites[jumps_start_id-1].set_site_atom(sites[jumps_end_id-1].get_site_atom());
			sites[jumps_end_id-1].set_site_atom(atom_tmp);
			
			//startとendのdiffusion_idを交換
			int diff_tmp = sites[jumps_start_id-1].get_diffusion_id();
			sites[jumps_start_id-1].set_diffusion_id(sites[jumps_end_id-1].get_diffusion_id());
			sites[jumps_end_id-1].set_diffusion_id(diff_tmp);
			
			//diffusion_idをもつdiffusion_species[x].jump_totalにjump_vectorを追加する
			for (int i = 0; i != diffusion_species.size(); i++) {
				if ( diffusion_species[i].get_diffusion_id() == diff_tmp ) {
					vector<double> jump_total_tmp(3,0.0);
					for (int k = 0; k != jump_total_tmp.size(); k++) {
						jump_total_tmp[k] = diffusion_species[i].get_jump_total()[k] + jumps_possible[jump_happen_number].get_jump_vector()[k];
						
					}
					diffusion_species[i].set_jump_total(jump_total_tmp);


/*					//確認用
					for (int j = 0; j != diffusion_species[i].get_jump_total().size();j++)
					cout << "diffusion_species[" << i << "].jump_total = " << diffusion_species[i].get_jump_total()[j] << endl;
*/

				}

			}






			//時間を更新する
			double rho_2 = random0to1(mt);
			//cout << "\t" << "\t" << "rho_2 = " << rho_2 << endl;

			double dt = -log(rho_2)/freq_sum;
			total_time += dt;
			//cout << "dt = " << dt << endl;
			//cout << "\t" << "\t" << "total_time = " << total_time << endl;
			//cout << endl;

/*			//確認用
			for (int i = 0; i != sites.size(); i++) {
				cout << "site[" << i << "] ";
				cout << "site_id = " << sites[i].get_site_id() << " " ;
				cout << "site_atom = " << sites[i].get_site_atom() << " " ;
				cout << "diffusion_id = " << sites[i].get_diffusion_id() << endl;
			}
*/


			//実行時間の計測
	//		end = chrono::system_clock::now();

	//		auto time = end - start;

	//		time_stamp = chrono::system_clock::to_time_t(start);
			//cout << "\t" << ctime(&time_stamp);

	//		auto msec = chrono::duration_cast<chrono::microseconds>(time).count();
			//cout << "\t" << msec << " msec" << endl;
			




		}

		//ループ終了後
		//cout << "\t" << "loop_counter finished" << endl;
		//cout << endl;
		//
		//拡散係数を出力およびファイルに出力する


		//自己拡散係数を求めるために拡散種ごとのjump_totalを準備する
		vector<double> jump_total_all(3,0.0);

		//拡散種1つ1つに対し操作を行う
		for (int j = 0; j != diffusion_species.size(); j++) {


			//自己拡散係数(および電場勾配化では伝導度)を求めるため
			//拡散種ごとのjump_totalを足して合計変位を出す
			for (int i = 0; i != jump_total_all.size(); i++) {
				jump_total_all[i] += diffusion_species[j].get_jump_total()[i];
			}

			if (!E_field_yes) {   //電場なしのとき 

			//トレーサー拡散係数を計算し出力
			vector<double> D_t_3d = diffusion_species[j].get_D(lattice_matrix,total_time);
			D_t_3d_vector.push_back(D_t_3d);

			//もともとOUTPUTに出力していたが不必要になった
			/*
			for (int i = 0; i != D_t_3d.size(); i++) {
				//cout << "\t"  << "diffusion_species[" << j << "]:D_t_3d[" << i << "] = " << D_t_3d[i] << endl;
				ofs << D_t_3d[i] << "," ;
				if (i+1 == D_t_3d.size()) {
					ofs << endl;
				}
			}
			//cout << endl;
			*/

			}



		} 

			
		//jump_total_allを分率座標からcartesian座標に直す(関数により用済みとなった)
		/*
		Eigen::Vector3d jump_total_all_cartesian;
		jump_total_all_cartesian << jump_total_all[0], jump_total_all[1], jump_total_all[2];
		Eigen::Vector3d displacement_vector = lattice_matrix*jump_total_all_cartesian;
		cout << "displacement_vector = " << displacement_vector << endl;
		*/

		//transcoords関数でjump_total_allを分率座標からcartesian座標に直す
		Eigen::Vector3d displacement_vector;
		displacement_vector = transcoords(jump_total_all,lattice_matrix);


		//電場ありの場合、伝導度テンソルを出力
		if (E_field_yes) {
			vector<double> Sigma_x(3,0.0);
			double concentration = diffusion_species.size() / lattice_matrix.determinant();
			//cout << "lattice_matrix.determinant() = " << lattice_matrix.determinant() << endl;
			//cout << "concentration = " << concentration << endl;
		
			for (int j = 0; j != Sigma_x.size(); j++) {
				
				Sigma_x[j] = displacement_vector[j];
				Sigma_x[j] *= q_charge * concentration / (E_field_strength * total_time);
				Sigma_x[j] *= pow(10,8); //Å^-1をcm^-1に変換
				//cout << "Sigma_x[" << j << "]  [S/cm] = " << Sigma_x[j] << endl;


			}
			
			Sigma_x_vector.push_back(Sigma_x);


		}



		//電場なしの場合、自己拡散係数を出力
		else {
			displacement_vector *= pow(10,-8);

			//jump_total_allを2乗したあと拡散種の数で割って、自己拡散係数を出力する
			vector<double> D_j_3d(3,0.0);

			for (int j = 0; j != jump_total_all.size(); j++) {
				
				D_j_3d[j] = pow(displacement_vector(j), 2.0);
				D_j_3d[j] /= diffusion_species.size();
				D_j_3d[j] /= 2*total_time;
				//cout << "\t"  << "D_j_3d[" << j << "] = " << D_j_3d[j] << endl;
			}

			D_j_3d_vector.push_back(D_j_3d);
		}



		//cout << endl;

		//sitesをリセットする
		for (int i = 0; i != sites.size(); i++ ) {
			sites[i].set_site_atom(1);
			sites[i].set_diffusion_id(-1);
		}

		//実行時間の計測
	//	end = chrono::system_clock::now();

	//	auto time = end - start;

	//	time_stamp = chrono::system_clock::to_time_t(start);
		//cout <<  ctime(&time_stamp);

	//	auto msec = chrono::duration_cast<chrono::microseconds>(time).count();
		//cout <<  msec << " msec" << endl;
		//cout << endl;

		//KMCが何回終わったか
		cout << '\t' <<  step_counter << " times KMC finished" << endl;
		ofs_log << '\t' <<  step_counter << " times KMC finished" << endl;


	}



	//電場なしの場合
	if (!E_field_yes) {

		//結果を出力するアウトプットファイルを作成する
		ofstream ofs_diff("DiffusionCoefficient", ios::app);

		//トレーサー拡散係数
		ofs_diff << "トレーサー拡散係数 (cm^2/s)" << endl;
		double Dx = 0;
		double Dy = 0;
		double Dz = 0;
		for (int i = 0; i != D_t_3d_vector.size(); i++) {

			Dx += D_t_3d_vector[i][0];
			Dy += D_t_3d_vector[i][1];
			Dz += D_t_3d_vector[i][2];
				
			if (i+1 == D_t_3d_vector.size()) {

				Dx /= D_t_3d_vector.size();
				Dy /= D_t_3d_vector.size();
				Dz /= D_t_3d_vector.size();
			}

		}

		ofs_diff << "Dx = " << Dx << endl; 
		ofs_diff << "Dy = " << Dy << endl; 
		ofs_diff << "Dz = " << Dz << endl; 
		
		
		//自己拡散係数
		ofs_diff << "自己拡散係数 (cm^2/s)" << endl;
		Dx = 0;
		Dy = 0;
		Dz = 0;
		for (int i = 0; i != D_j_3d_vector.size(); i++) {

			Dx += D_j_3d_vector[i][0];
			Dy += D_j_3d_vector[i][1];
			Dz += D_j_3d_vector[i][2];
				
			if (i+1 == D_j_3d_vector.size()) {

				Dx /= D_j_3d_vector.size();
				Dy /= D_j_3d_vector.size();
				Dz /= D_j_3d_vector.size();
			}

		}

		ofs_diff << "Dx = " << Dx << endl; 
		ofs_diff << "Dy = " << Dy << endl; 
		ofs_diff << "Dz = " << Dz << endl; 

		//化学拡散係数
		ofs_diff << "化学拡散係数 (cm^2/s)" << endl;

	}



	//電場ありの場合
	if (E_field_yes) {

		//結果を出力するアウトプットファイルを作成する
		ofstream ofs_sigma("ElectricalConductivity", ios::app);
		ofstream ofs_sigma_2("sigma_log", ios::app);
		

		//伝導度x成分
		ofs_sigma << "伝導度 [S/cm]" << endl;
		vector<double> Sigma_x_vector_total(3,0.0);
		for (int j = 0; j != Sigma_x_vector_total.size(); j++) {
			for (int i = 0; i != Sigma_x_vector.size(); i++) {

				Sigma_x_vector_total[j] += Sigma_x_vector[i][j];
				ofs_sigma_2 << "Sigma_x_vector[" << i << "][" << j << "] = " << Sigma_x_vector[i][j] << endl;

				if (i+1 == Sigma_x_vector.size()) {
					Sigma_x_vector_total[j] /= Sigma_x_vector.size();
					ofs_sigma << "Sigma_x" << j << " = " << scientific << Sigma_x_vector_total[j] << " ,  " ;
					ofs_sigma_2 << endl;

				}
			
			}

			double variance_sigma = 0.0 ;
			double standard_deviation = 0.0 ;
			for (int i = 0; i != Sigma_x_vector.size(); i++) {
				variance_sigma += pow(Sigma_x_vector[i][j] - Sigma_x_vector_total[j], 2);

				if (i+1 == Sigma_x_vector.size()) {
					variance_sigma /= Sigma_x_vector.size();
					standard_deviation = sqrt(variance_sigma);
					ofs_sigma << "Standard deviation = " << scientific << standard_deviation << endl;
					
				}

			}
		}

		

		
		



	}

	//実行時間の計測
	end = chrono::system_clock::now();
	auto time = end - start;

	time_stamp = chrono::system_clock::to_time_t(start);
	ofs_log <<  endl;

	auto sec = chrono::duration_cast<chrono::seconds>(time).count();

	ofs_log << "Execution time = " << sec << " sec" << endl;
	cout << endl;
	
	
}

