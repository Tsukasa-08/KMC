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
#include <set>
#include <map>
#include <numeric>
#include <ctime>
#include <iomanip>

#include "Site.h"
//#include "Jump.h"
#include "Diffusionspecie.h"
#include "param.hpp"
#include "toml.hpp"

using namespace std;

//
/*
//サイトクラスを定義
class Site {

private:
	int site_id;
	int site_atom;
	int diffusion_id;
	std::vector<double> site_frac_coords;
	std::vector<Jump> jumps_from_here;
	std::vector<Jump> jumps_to_here;
	std::vector<int> blocking_mate_list;

public:
	
	//デフォルトコンストラクタ デフォルトでは空孔にしておく
	Site() : site_id(0), site_frac_coords(3,0.0), site_atom(1), diffusion_id(-1){
	}

	//コピーコンストラクタ
	Site(const Site &src){
	}

	//セッタ
	void set_site_id(int id) { site_id = id; };

	void set_site_atom(int atom) { site_atom = atom; };

	void set_diffusion_id(int d_id) { diffusion_id = d_id; };

	void set_site_frac_coords(std::vector<double> v) { site_frac_coords = v; };

	void set_jumps_from_here(std::vector<Jump> here) { jumps_from_here = here ; } ;

	void set_a_jump_jumps_from_here(Jump here) { jumps_from_here.push_back(here) ; } ;

	void set_a_jump_jumps_to_here(Jump here) { jumps_to_here.push_back(here) ; } ;

	void set_blocking_mate_list(std::vector<int> st) { blocking_mate_list = st ; } ;

	//ゲッタ
	int get_site_id() { return site_id; } ;

	int get_site_atom() { return site_atom; } ;

	int get_diffusion_id() { return diffusion_id; } ;

	std::vector<double> get_site_frac_coords() { return site_frac_coords; } ;

	std::vector<Jump> get_jumps_from_here() { return jumps_from_here; } ;

	std::vector<Jump> get_jumps_to_here() { return jumps_to_here; } ;

	std::vector<int> get_blocking_mate_list() { return blocking_mate_list; } ;
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
	int diffusion_siteid_now;
	std::vector<double> jump_total;
	std::vector<double> sum_squared_distance;

public:
	int diffusion_counter;
	static std::vector<int> diffusion_siteid_now_list;
	static std::set<int> blocking_list;

	Diffusionspecie() : diffusion_id(-1), jump_total(3,0.0) {
	};

	//セッタ
	void set_diffusion_id(int d_id) { diffusion_id = d_id; }

	void set_diffusion_siteid_now(int d_sid_now) { diffusion_siteid_now = d_sid_now; }

	void set_diffusion_counter(int d_ctr) { diffusion_counter = d_ctr; }

	void set_jump_total(std::vector<double> jt) { jump_total = jt; }

	void set_sum_squared_distance(std::vector<double> sum_s_d) { sum_squared_distance = sum_s_d; }

	void set_diffusion_siteid_now_list(std::vector<int> ds_list) { diffusion_siteid_now_list = ds_list; }

	//ゲッタ
	int get_diffusion_id() { return diffusion_id; }

	int get_diffusion_siteid_now() { return diffusion_siteid_now; }

	int get_diffusion_counter() { return diffusion_counter; }

	std::vector<double> get_jump_total() { return jump_total; }

	std::vector<double> get_sum_squared_distance() { return sum_squared_distance; }

	std::vector<int> get_diffusion_siteid_now_list() { return diffusion_siteid_now_list; }

	
	//拡散係数を算出する関数(引数にlattice_matrixとtotal_time)
	std::vector<double> get_D(Eigen::Matrix3d lattice_matrix, double t) {
		Eigen::Vector3d eigen_jump_total;
		eigen_jump_total << get_jump_total()[0], get_jump_total()[1], get_jump_total()[2];
		Eigen::Vector3d displacement_vector = lattice_matrix*eigen_jump_total;
		displacement_vector *= pow(10,-8);
		std::vector<double> D_3d(3,0.0);
		for (int i = 0; i != D_3d.size(); i++) { 
			double mean_square_displacement = pow(displacement_vector(i),2.0);
			D_3d[i] = mean_square_displacement/(2*t);
		}
		return D_3d;
	}

	//二乗変位を算出する関数(引数にlattice_matrix)
	std::vector<double> get_mean_square_displacement(Eigen::Matrix3d lattice_matrix) {
		Eigen::Vector3d eigen_jump_total;
		eigen_jump_total << get_jump_total()[0], get_jump_total()[1], get_jump_total()[2];
		Eigen::Vector3d displacement_vector = lattice_matrix*eigen_jump_total;
		displacement_vector *= pow(10,-8);
		std::vector<double> mean_square_displacement(3,0.0);
		for (int i = 0; i != mean_square_displacement.size(); i++) { 
			mean_square_displacement[i] = pow(displacement_vector(i),2.0);
		}
		return mean_square_displacement;
	}

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

//Eigen::VectorXdをstd::vectorに型変換
vector<double> eigen2vector(Eigen::Vector3d eigen_vector) 
{
	vector<double> vector_cartesian(3,0.0);
	for (int i = 0 ; i != vector_cartesian.size() ; i++) {
		vector_cartesian[i] = eigen_vector(i) ; 
	}

	return vector_cartesian;
}

//ある値がvector内の要素に含まれているか否か判定する
int vector_finder(std::vector<int> vec, int number) {
	auto itr = std::find(vec.begin(), vec.end(), number);
	size_t index = std::distance(vec.begin(), itr);
	if (index != vec.size()){
		return 1;
	}
	else {
		return 0;
	}
}




//原子種の数字を定義

int const number_vacancy = 1;
int const number_proton = 2;


//物理定数を定義

double const kb = 1.380649 * pow(10,-23);
double const ec = 1.60217662 * pow(10,-19);

//拡散種の電荷を定義(今回はプロトン=1),要改善
double ion_charge = 1;
double q_charge = ion_charge * ec;

//Nernst-Einsteinの関係式より、拡散係数から伝導度を算出する関数(1次元、1方向)
double NernstEinstein_DtoSigma(double D, double concentration, int temperture, double ion_charge, int dimensionality) {
	double Sigma; //[S/cm]
	Sigma = D * pow(ion_charge * ec, 2) * concentration / (kb * temperture) * pow(10, 8*dimensionality);
	return Sigma;

}

//Nernst-Einsteinの関係式より、伝導度から拡散係数を算出する関数(1次元、1方向)
double NernstEinstein_SigmatoD(double Sigma, double concentration, int temperture, double ion_charge, int dimensionality) {
	double D; //[cm^2/s]
	D = Sigma * (kb * temperture) * pow(10, -8*dimensionality) /  (pow(ion_charge * ec, 2) * concentration) ;
	return D;

}


//拡散種の配置一覧vectorを定義
std::vector<int> Diffusionspecie::diffusion_siteid_now_list;
std::set<int> Diffusionspecie::blocking_list;

//メイン関数の開始
int main()
{
	//プログラム開始時刻を表示
	chrono::system_clock::time_point start, end;
	chrono::milliseconds start_msec, end_msec;
	time_t time_stamp;
	start = chrono::system_clock::now();
	
	//開始時間のmsecを求める
	start_msec = chrono::duration_cast<chrono::milliseconds>(start.time_since_epoch());
	long long all_msec = start_msec.count();
	int msec = all_msec % 1000 ;

	time_stamp = chrono::system_clock::to_time_t(start);
	struct tm* timer = localtime(&time_stamp);
	cout << endl;
	cout << "#################################" << endl;
	cout << timer->tm_year + 1900 << "-" 
		<< timer->tm_mon +1 << "-" 
		<< timer->tm_mday << " " 
		<< timer->tm_hour << ":"
		<< timer->tm_min << ":"
		<< setfill('0') << right << setw(2) << timer->tm_sec << setfill(' ') << "." 
		<< setfill('0') << right << setw(3) << msec << setfill(' ') << endl;
	cout << "Start program." << endl;
	cout << "#################################" << endl;
	cout << endl;

	//アウトプットファイルを作成しておく
	ofstream ofs_output("OUTPUT", ios::app);

	//アウトプットファイル用のディクショナリやvectorを作成しておく
	map<string, double> map_for_output;
	vector<double> total_time_for_output;

	//変位一覧を出力するmean_displacement.csvを開いておく
	ofstream ofs_ave_dis("mean_displacement.csv", ios::app);
	//ofs_ave_dis << "#KMC " << step_counter << " times" << endl;
	ofs_ave_dis << "#the number of KMC, diffusion_id, mean displacement in x [Ang.], mean_displacement in y [Ang.], mean displacement in z [Ang.], sum of squared displacement of each jumps in x, y, z [Ang.^2], jump_counter [times]" << endl;
	ofs_ave_dis << "KMC_times,diffusion_id,dx,dy,dz,sum_x2,sum_y2,sum_z2,counter" << endl;

	//結果を出力するアウトプットファイルを作成する
	ofstream ofs_diff("DiffusionCoefficient", ios::app);

	//結果を出力するアウトプットファイルを作成する
	ofstream ofs_sigma("IonicConductivity", ios::app);

	//実行時間の計測

	//make DiffusionCoefficient vector
	vector< vector<double> > D_t_3d_vector;
	vector< vector<double> > D_j_3d_vector;
	vector< vector<double> > D_c_3d_vector;

	//make ElectricalConductivity vector
	vector< vector<double> > Sigma_vector;

	//make average_displacement vector
	vector< vector<double> > average_displacement_vector;

	//INPUTを読み込む
	param::parameter param("INPUT");
	long long mcsp = param.get<int>("MCSP", 0);
	int average = param.get<int>("AVERAGE", 0);
	long long p_place_n = param.get<int>("NDIFFS", 0);
	int E_field_yes = param.get<int>("EFIELDON", 0);
	//double E_field_strength_for_pow = param.get<double>("EFIELD", 0);
	double correct_constant_for_pow = param.get<double>("CORRECT", 0);
	double temperture = param.get<double>("TEMP", 0);
	int E_field_axis = param.get<int>("AXIS", 0);
	double distance_jump = param.get<double>("DISTANCEJUMP", 1); //単位は[Å]
	int blocking_yes = param.get<int>("BLOCKING", 0);
	int dimensionality = param.get<int>("DIMENSIONALITY", 3);

	//読み込んだINPUTをもとに計算
	int step_max = (average + p_place_n - 1) / p_place_n;
	long long loop_max = mcsp * p_place_n;
	double correct_constant = pow(10, correct_constant_for_pow);
	//double E_field_strength = pow(10, E_field_strength_for_pow);
	double E_field_strength = correct_constant * kb * temperture / (q_charge * 0.5 * distance_jump); //単位は[J/Å]
	
	//次元に応じて伝導度の単位を対応させておく
	string conductivity_unit;
	switch (dimensionality) {
		case 1 : conductivity_unit = "[S*cm]"; break;
		case 2 : conductivity_unit = "[S]"; break;
		case 3 : conductivity_unit = "[S/cm]"; break;
	}
	
	//後半濃度算出用
	double concentration;


	
	//入力確認用logファイルを作成

	//読み込めたか確認用
	if (!param) {
		cerr << "Could not find file INPUT" << endl;
		abort();
	}
	cout << "INPUT read" << endl;
	cout << "\t" << "Monte Carlo Step per Particle : MCSP = " << mcsp << endl;
	cout << "\t" << "the number of particles for calculate ensemble average :  AVERAGE = " << average << endl;
	cout << "\t" << "the number of executing KMCs : NSTEPS = " << step_max << endl;
	cout << "\t" << "Monte Carlo Step (MCS) :  NLOOPS = " << loop_max << endl;
	cout << "\t" << "the number of placing diffusion atoms :  NDIFFS = " << p_place_n << endl;
	cout << "\t" << "Electrical field strength :  EFIELD = " << scientific << E_field_strength << " [V/Å] " << endl;
	cout << "\t" << "Electrical field strength :  EFIELD = " << E_field_strength*pow(10,8) << defaultfloat << " [V/cm] " << endl;
	cout << "\t" << "temperture : TEMP = " << temperture << endl;
	cout << "\t" << "the correct_constant for adjusting jump frequency gamma (default = 0.1) :  correct_constant = " << correct_constant << endl;
	cout << "\t" << "the direction of E_field : E_field_axis = " << E_field_axis << " (+x=1, +y=2, +z=3, -x=-1, -y=-2, -z=-3)" << endl;
	cout << "\t" << "whether blocking_list is valid or not = " << blocking_yes << " (1=valid、0=invalid)"<< endl;
	cout << "\t" << "the dimensionality in this diffusion situation : dimensionality = " << dimensionality << endl;
	cout << endl;

	
	//電場の向きにより出力する伝導度を変える
	string axis;
	switch (E_field_axis) {
		case 1 : axis = "x"; break;
		case 2 : axis = "y"; break;
		case 3 : axis = "z"; break;
		case -1 : axis = "-x"; break;
		case -2 : axis = "-y"; break;
		case -3 : axis = "-z"; break;
	}

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
		}
	}

	//サイト数=POSCARの行数-(座標一覧までの行数)
	site_total_number -= DIRECT_num;
	cout << "POSCAR read" << endl;
	cout << "\t" << "string [Direct] is on the " << DIRECT_num << " line." << endl;
	cout << "\t" << "site_total_number = " << site_total_number << endl;
	cout << endl;


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
			for (int i = 0; i <= 2; i++) {
				ss_line >> s;
				frac_dblvec[i] = stod(s);
				//cout << "frac_dblvec[" << i << "] = "<< frac_dblvec[i] << endl;
				
			}

			//1行を空白で3つの座標に分割していく
			//while(getline(ss_line, s, ' ')){

			//	frac_dblvec[s_counter] = stod(s);
			//	cout << "frac_dblvec[" << s_counter << "] = "<< frac_dblvec[s_counter] << endl;
			//	s_counter += 1;
			//}

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


/*	//Site確認用
	for (int i = 0; i != sites.size() ; i++) {
		cout << "sites[" << i << "].id = " << sites[i].get_site_id() << endl;
		for (int j = 0; j <= 2; j++){
			cout << '\t' << sites[i].get_site_frac_coords()[j] << endl;
		}
	}

	
*/



	//blocking_list.csvを読み込む
	int csv_total_number = 0;
	if (blocking_yes) {
		ifstream for_line3("blocking_list.csv");

		//読み込めたか確認用
		if (!for_line3) {
			cerr << "Could not find file blocking_list.csv" << endl;
			abort();
		}

		string for_line_reader3;
		while (getline(for_line3,for_line_reader3)) {

			csv_total_number++;
		}

		cout << "blocking_list.csv read" << endl;
		cout << endl;
	}

	vector< vector<int> > blocking_list_csv(csv_total_number);

	if (blocking_yes) {
	ifstream ifs3("blocking_list.csv");
	n_lines = 0;
	while (getline(ifs3, line)){


		vector<string> strvec = split(line, ',');

		for (int i = 0; i != strvec.size(); i++){
			blocking_list_csv[n_lines].push_back(stoi(strvec[i]));
		}


	//	cout << "jumps[" << n_lines << "].freq = " << jumps[n_lines].get_freq() << endl;
		
		n_lines += 1;
	}	

/*	//csv確認用
	for (auto itr = blocking_list_csv.begin(); itr != blocking_list_csv.end(); itr++) {
		for (auto itr2 = (*(itr)).begin(); itr2 != (*(itr)).end(); itr2++){
			cout << *itr2 << "," ;
		}
		cout << endl;
	}

*/

	
	//blocking_list_csvをもとに、各Siteのblocking_mate_listに追加していく
	for (auto itr = blocking_list_csv.begin(); itr != blocking_list_csv.end(); itr++) {
		for (auto itr2 = (*(itr)).begin(); itr2 != (*(itr)).end(); itr2++){
			sites[*itr2-1].set_blocking_mate_list(*itr);
		}
	}

	/*//mateがSitesに属しているか確認用
	for (int i = 0; i != sites.size(); i++) {
		cout << "sites " << sites[i].get_site_id() << "のmate一覧" << endl;
		vector<int> mate_tmp = sites[i].get_blocking_mate_list();
	
		for (auto itr = mate_tmp.begin() ; itr != mate_tmp.end() ; itr++) {
			cout << '\t' << *itr << endl;
		}	

	}
	*/
	

	}
		
	



	//Site.site_frac_coordsをもとに、Jump.jump_vectorを生成する
	for (int i = 0; i != jumps.size(); i++) {

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
		}
		

		//jump_vector_tmpをjump_vectorに代入する
		jumps[i].set_jump_vector(jump_vector_tmp);


/*		//確認用
		for (int j = 0; j != jumps[i].get_jump_vector().size(); j++) {
			cout << "jump[" << i << "].jump_vector[" << j << "] = " << jumps[i].get_jump_vector()[j] << endl;
		}
		*/
		

		


	}

/*//ジャンプがSitesに属しているか確認用
	for (int i = 0; i != sites.size(); i++) {
		cout << "sites[" << i << "]のjump一覧" << endl;
	
		for (int j = 0; j != sites[i].get_jumps_to_here().size() ; j++) {
			cout << "vector " << j << "番目" << endl;
		

			for (int k = 0 ; k <= 2 ; k++) {

				cout << '\t' << sites[i].get_jumps_to_here()[j].get_jump_vector()[k] << endl;
			}

		}	

	}
	*/




	//電場がかかっていた場合、ジャンプ頻度を補正する
	if (E_field_yes) {

		//cartesian座標軸方向の単位ベクトルを作成し、電場の大きさをかけて電場ベクトルとする
		//変数E_field_axisによって作成する単位ベクトルを変える(x=0, y=1, z=2)
		
		Eigen::Vector3d E_field_vector(3); 

		switch (E_field_axis) {
			case 1 : E_field_vector << 1,0,0 ; break;
			case 2 : E_field_vector << 0,1,0 ; break;
			case 3 : E_field_vector << 0,0,1 ; break;
			case -1 : E_field_vector << -1,0,0 ; break;
			case -2 : E_field_vector << 0,-1,0 ; break;
			case -3 : E_field_vector << 0,0,-1 ; break;
		}


		E_field_vector *= E_field_strength;
		cout << "E_field_vector = [ " << E_field_vector(0) << ", " << E_field_vector(1) << ", " << E_field_vector(2)  << " ]" << endl;

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
		//	cout << "jump_frep = " << jumps[i].get_freq() << endl;
		//	cout << "fixed_jump_freq = " << fixed_jump_freq << endl;
			jumps[i].set_freq(fixed_jump_freq);
		}

	}

	else {
		cout << "no E_field" << endl;	
	}

	//補正したジャンプ頻度情報をSiteと結びつける
	for (int i = 0, n = jumps.size(); i != n; i++) {
		
		//始点と終点のsite_idを取得
		int start_site_id_tmp = jumps[i].get_start_site_id();
		int end_site_id_tmp = jumps[i].get_end_site_id();
		
		//jump_vector_tmpをSites.jumps_from_hereにpush_backして格納する(set_以下はpush_back用の関数)
		sites[start_site_id_tmp-1].set_a_jump_jumps_from_here(jumps[i]);

		sites[end_site_id_tmp-1].set_a_jump_jumps_to_here(jumps[i]);
	}


/*	//確認用
	for (int i = 0; i != sites.size(); i++) {
		for (int j = 0; j != sites[i].get_jumps_from_here().size(); j++) {

			cout << "\t" << "\t" << "start = " << sites[i].get_jumps_from_here()[j].get_start_site_id() ; 
			cout << "\t" << "\t" << "end = " << sites[i].get_jumps_from_here()[j].get_end_site_id() ; 
			cout << "\t" << "\t" << "jumps_from_here[" << j << "] = " << sites[i].get_jumps_from_here()[j].get_freq() << endl;
		}
	}
	cout << endl;
	//
	//確認用
	for (int i = 0; i != sites.size(); i++) {
		for (int j = 0; j != sites[i].get_jumps_to_here().size(); j++) {

			cout << "\t" << "\t" << "start = " << sites[i].get_jumps_to_here()[j].get_start_site_id() ; 
			cout << "\t" << "\t" << "end = " << sites[i].get_jumps_to_here()[j].get_end_site_id() ; 
			cout << "\t" << "\t" << "jumps_to_here[" << j << "] = " << sites[i].get_jumps_to_here()[j].get_freq() << endl;
		}
	}

	return 0;
	*/


	concentration = p_place_n / lattice_matrix.determinant();
	


	for (int step_counter = 1; step_counter <= step_max; step_counter++) {
		
		//KMC開始時刻を表示
		chrono::system_clock::time_point kmc_start, kmc_end;
		chrono::milliseconds kmc_start_msec, kmc_end_msec;
		time_t kmc_time_stamp;
		kmc_start = chrono::system_clock::now();
		
		//開始時間のmsecを求める
		kmc_start_msec = chrono::duration_cast<chrono::milliseconds>(kmc_start.time_since_epoch());
		long long kmc_all_msec = kmc_start_msec.count();
		int kmc_msec = kmc_all_msec % 1000 ;

		time_stamp = chrono::system_clock::to_time_t(kmc_start);
		struct tm* kmc_timer = localtime(&time_stamp);
		cout << endl;
		cout << "\t" << "#################################" << endl;
		cout << "\t" << kmc_timer->tm_year + 1900 << "-" 
			<< kmc_timer->tm_mon +1 << "-" 
			<< kmc_timer->tm_mday << " " 
			<< kmc_timer->tm_hour << ":"
			<< kmc_timer->tm_min << ":"
			<< setfill('0') << right << setw(2) << kmc_timer->tm_sec << setfill(' ') << "." 
			<< setfill('0') << right << setw(3) << kmc_msec << setfill(' ') << endl;
		cout << "\t" << "Start KMC " << step_counter << " times" << endl;
		cout << "\t" << "#################################" << endl;
		cout << endl;





		//確認用
		//cout << "step_counter = " << step_counter << " start " << endl;

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
		vector<int> random_proton_place_number_vector;

		//blockingありの場合
		if (blocking_yes) {
			//cout << "blockingあり" << endl;

			//乱数を生成し、blocking_list_csv(vector<vector<int>>)をシャッフルする
			shuffle( blocking_list_csv.begin(), blocking_list_csv.end(), mt );
			
			//先頭から抽出する
			for (int i = 0; i != p_place_n; i++) {
				shuffle ( blocking_list_csv[i].begin(), blocking_list_csv[i].end(), mt);
				random_proton_place_number_vector.push_back(blocking_list_csv[i][0]);
			}

		}

		//blockingなしの場合
		else {
			//cout << "blockingなし" << endl;
			//(サイトの数-1)を要素にもつvectorを生成
			random_proton_place_number_vector.resize(site_total_number-1);
			for (int i = 0; i != random_proton_place_number_vector.size(); i++) {
				random_proton_place_number_vector[i] = i + 1;
			}

			//乱数を生成し、vectorをシャッフルする
			shuffle( random_proton_place_number_vector.begin(), random_proton_place_number_vector.end(), mt );
			

			//先頭からプロトンを配置する数分(p_place_n)だけ抜き出しソートする
			random_proton_place_number_vector.resize(p_place_n);
			sort( random_proton_place_number_vector.begin(), random_proton_place_number_vector.end() );

		}

		/*//確認用
		for (int i = 0; i != random_proton_place_number_vector.size(); i++) {
			cout << "random_proton_place_number_vector[" << i << "] = " << random_proton_place_number_vector[i] << endl;
		}

		return 0;
		*/
		
		

		//Diffusionspecieクラスのvectorをつくる(数はプロトン配置数=p_place_n)
		vector<Diffusionspecie> diffusion_species(p_place_n);
		//それぞれのidを[i]にたいしてi+1で設定する(1からp_place_nまでdiffusion_idとして通し番号をふる)
		for (int i = 0; i != diffusion_species.size(); i++) {
			diffusion_species[i].set_diffusion_id(i+1);
			diffusion_species[i].diffusion_counter = 0;
		}
		
		


		//プロトンを配置するSiteクラスのsite_atomをnumber_proton(今回は2を割当)に変更する
		//と同時に、diffusion_idを設定する
		
		int d_id = 1;

		for (int i = 0; i != random_proton_place_number_vector.size(); i++) {
			for (int k = 0; k != sites.size(); k++) {

				//sites[k]=Siteクラスのsite_idが、プロトンを配置するsite_idかどうかを判定
				if (sites[k].get_site_id() == random_proton_place_number_vector[i]) {
					sites[k].set_site_atom(number_proton);
					sites[k].set_diffusion_id(d_id);

					//同時に、静的メンバ変数であるdiffusion_siteid_now_listにも追加しておく
					Diffusionspecie::diffusion_siteid_now_list.push_back(sites[k].get_site_id());

					//同時に、静的メンバ変数であるblocking_listにも追加しておく
					vector<int> temp_for_add_blocking_list = sites[k].get_blocking_mate_list();
					for (auto itr = temp_for_add_blocking_list.begin() ; itr != temp_for_add_blocking_list.end() ; itr++) {
						Diffusionspecie::blocking_list.insert(*itr);
					}


					d_id++;

					//cout << "proton at site id = " << sites[k].get_site_id() << endl;
					break;
				}
				else {
					continue;
				}
			}
		}

/*		//確認用
		cout << "blocking_list = " ;
		for (auto itr = Diffusionspecie::blocking_list.begin(); itr != Diffusionspecie::blocking_list.end() ; itr++) {
			cout << *itr << "," ;
		}
		cout << endl;
*/





	//	if ( d_id != p_place_n+1 )
			//cout << "d_id != p_place_n+1" << endl;
	//	else
			//cout << "d_id = p_place_n+1" << endl;

		//確認用
/*		for (int i = 0; i != sites.size(); i++) {
			cout << "site[" << i << "] ";
			cout << "site_id = " << sites[i].get_site_id() << " " ;
			cout << "site_atom = " << sites[i].get_site_atom() << " " ;
			cout << "diffusion_id = " << sites[i].get_diffusion_id() << endl;
		}
		*/

		
		
		//シミュレーション時間total_timeを定義
		double total_time = 0.0;


		//メインループの開始
		cout << endl;
		
		int start_was = 0;
		int end_was = 0;
		vector<Jump> jumps_possible_tmp;
		vector<Jump> jumps_possible;
		vector<Jump> jumps_impossible_tmp;

		//cout << "\t" << "main loop start" << endl;
		//cout << "\t" << "loop processing…" << endl;
		string processing;
		for (long long loop_counter = 1; loop_counter <= loop_max; loop_counter++) { 			


			//実行時間を計測する
	//		chrono::system_clock::time_point start, end;
	//		time_t time_stamp;

	//		start = chrono::system_clock::now();


			//確認用
			//cout << endl;
			//cout << "\t" << "loop_counter = " << loop_counter << endl;
			
			if (loop_counter % (loop_max/10) == 0) {
				processing += "#";
				cout  << "\t" << processing << endl;
				if (loop_counter == loop_max ) {
					cout << endl;
				}
			}
			



			//1ループ目のみ
			if (loop_counter == 1) {
						
				//系で起きうる事象(今回はジャンプ)を列挙し、jumps_possibleに入れていく
				//vector<Jump> jumps_possible_tmp;
				//vector<Jump> jumps_possible;

				//diffusion_siteid_now_listにあるサイト上からのジャンプを追加したい
				for (int i = 0, n = Diffusionspecie::diffusion_siteid_now_list.size() ; i != n; i++) {
					
					//そのサイトからのジャンプを取得
					vector<Jump> temp_jump_vector_loop_1 = sites[Diffusionspecie::diffusion_siteid_now_list[i]-1].get_jumps_from_here();

					//blockingあり
					if (blocking_yes) {
						for (auto itr = temp_jump_vector_loop_1.begin(); itr != temp_jump_vector_loop_1.end(); itr++) {

							//回転経路=同一blocking area内でのジャンプの場合は追加する
							if (vector_finder(sites[Diffusionspecie::diffusion_siteid_now_list[i]-1].get_blocking_mate_list(), (*itr).get_end_site_id())) {
								jumps_possible.push_back(*itr);
							}
							//ホッピング経路=別のblocking areaへのジャンプの場合はblocking_listにプロトンがいなければ追加
							else {
								if (!Diffusionspecie::blocking_list.count((*itr).get_end_site_id())) {
									jumps_possible.push_back(*itr);
								}
							}
						}
						
					}

					//blockingなし
					else {
						//終点にプロトンがいなければ追加
						for (auto itr = temp_jump_vector_loop_1.begin(); itr != temp_jump_vector_loop_1.end(); itr++) {
							if (!vector_finder(Diffusionspecie::diffusion_siteid_now_list, (*itr).get_end_site_id())) {
								jumps_possible.push_back(*itr);
							}
						}
					}

					
				}
				
				/*for (int i = 0, n = jumps.size() ; i != n; i++) {
					
					//jumps[i]の始点サイトidが、プロトンの現在サイトのリストに入っていればtmpに追加
					if (vector_finder(Diffusionspecie::diffusion_siteid_now_list, jumps[i].get_start_site_id())) {
						jumps_possible_tmp.push_back(jumps[i]);
					}
				}

				for (int i = 0, n = jumps_possible_tmp.size() ; i != n; i++) {
					
					//jumps[i]の終点サイトidが、プロトンの現在サイトのリストに入っていなければtmpに追加
					if (!vector_finder(Diffusionspecie::diffusion_siteid_now_list, jumps_possible_tmp[i].get_end_site_id())) {
						jumps_possible.push_back(jumps_possible_tmp[i]);
					}
					else {
						jumps_impossible_tmp.push_back(jumps_possible_tmp[i]);
					}
				}
				*/


			}


			//2ループ目以降
			else {

				//jumps_possibleを更新する

				//blockingあり
				if (blocking_yes) {
					
					auto itr = jumps_possible.begin();
					while (itr != jumps_possible.end()) {

						//start_wasからのjumpを削除
						if ((*itr).get_start_site_id() == start_was) {
							itr = jumps_possible.erase(itr);
						}

						//end_wasおよびそのblocking_mate_listに向かうjumpを削除

						else if (vector_finder(sites[end_was-1].get_blocking_mate_list(), (*itr).get_end_site_id())) {
							itr = jumps_possible.erase(itr);
						}

						//他はスルー
						else {
							itr++;
						}
					}

					//回転経路=同一blocking area内で移動した場合
					if (vector_finder(sites[start_was-1].get_blocking_mate_list(), end_was)) {

						//cout << "\t" << "rotation happend" << endl;

						/*//start_wasに向かうjumpを追加する…end_wasからのjumpですべてまかなえることに気づいた
						vector<Jump> temp_jump_vector_to_start = sites[start_was-1].get_jumps_to_here();
						for (auto itr_1 = temp_jump_vector_to_start.begin(); itr_1 != temp_jump_vector_to_start.end(); itr_1++) {

							//回転経路内かつ始点にプロトンがいれば追加
							if (vector_finder(sites[start_was-1].get_blocking_mate_list(), (*itr_1).get_start_site_id())
								&&
								vector_finder(Diffusionspecie::diffusion_siteid_now_list, (*itr_1).get_start_site_id())){
								jumps_possible.push_back(*itr_1);
							}
							//ホッピング経路なら追加はできないのでスルー
						}
						*/

						//end_wasからのjumpを追加する
						vector<Jump> temp_jump_vector_from_end = sites[end_was-1].get_jumps_from_here();
						for (auto itr_1 = temp_jump_vector_from_end.begin(); itr_1 != temp_jump_vector_from_end.end(); itr_1++) {
							
							//回転経路内なら追加
							if (vector_finder(sites[end_was-1].get_blocking_mate_list(), (*itr_1).get_end_site_id())) {
								jumps_possible.push_back(*itr_1);
							}

							//ホッピング経路の場合
							else {
								//終点がblocking_listになければ追加
								if (!Diffusionspecie::blocking_list.count((*itr_1).get_end_site_id())) {
									jumps_possible.push_back(*itr_1);
								}
							}
						}

					}

					//ホッピング経路=別のblocking areaへ移動した場合
					else {
						
						//cout << "\t" << "hopping happend" << endl;
						//cout << "\t" << "\t" << "search jump from end_was" << endl;

						//end_wasからのjumpを追加する
						vector<Jump> temp_jump_vector_from_end = sites[end_was-1].get_jumps_from_here();
						for (auto itr_2 = temp_jump_vector_from_end.begin(); itr_2 != temp_jump_vector_from_end.end(); itr_2++) {
							
							//回転経路内なら追加
							if (vector_finder(sites[end_was-1].get_blocking_mate_list(), (*itr_2).get_end_site_id())) {
								jumps_possible.push_back(*itr_2);
								//cout << "\t" << "\t" << "\t" << (*itr_2).get_start_site_id() << " to " << (*itr_2).get_end_site_id() << "  rot path add" << endl;
							}

							//ホッピング経路の場合
							else {
								//終点がblocking_listにない、かつ終点がstart_wasじゃなけば追加
								if (!Diffusionspecie::blocking_list.count((*itr_2).get_end_site_id())
									&&
									(*itr_2).get_end_site_id() != start_was) {
									//cout << "\t" << "\t" << "\t" << (*itr_2).get_start_site_id() << " to " << (*itr_2).get_end_site_id() << "  hop path add" << endl;
									jumps_possible.push_back(*itr_2);
								}
							}
						}

						//start_wasを含むblocking areaへのjumpを追加する
						//cout << "\t" << "\t" << "search jump to start_was" << endl;
						//blocking_areaのsite_id(vector)を取得
						vector<int> blocking_area_ids = sites[start_was-1].get_blocking_mate_list();
						//blocking_area内のsiteについて、それぞれに向かうjumpを探して追加する
						for (auto itr_mate_id = blocking_area_ids.begin(); itr_mate_id != blocking_area_ids.end(); itr_mate_id++) {
							vector<Jump> temp_jump_vector_to_start = sites[*itr_mate_id-1].get_jumps_to_here();
							for (auto itr_2 = temp_jump_vector_to_start.begin(); itr_2 != temp_jump_vector_to_start.end(); itr_2++){
								
								//回転経路内にプロトンは存在しないので何もしない
								if (vector_finder(sites[*itr_mate_id-1].get_blocking_mate_list(), (*itr_2).get_start_site_id())) {
								}

								//ホッピング経路の場合
								else {
									
									//始点にプロトンがあれば追加
									if (vector_finder(Diffusionspecie::diffusion_siteid_now_list, (*itr_2).get_start_site_id())) {
										//cout << "\t" << "\t" << "\t" << (*itr_2).get_start_site_id() << " to " << (*itr_2).get_end_site_id() << "  hop path add" << endl;
										jumps_possible.push_back(*itr_2);
									}
								}
							}
						}
					}
					
					
				}

				//blockingなし
				else {
					auto itr = jumps_possible.begin();
					while (itr != jumps_possible.end()) {

						//start_wasからのjumpを削除
						if ((*itr).get_start_site_id() == start_was) {
							itr = jumps_possible.erase(itr);
						}

						//end_wasに向かうjumpを削除
						else if ((*itr).get_end_site_id() == end_was) {
							itr = jumps_possible.erase(itr);
						}

						//他はスルー
						else {
							itr++;
						}
					}

					//start_wasに向かうjumpを追加
					vector<Jump> temp_jump_vector_to_start = sites[start_was-1].get_jumps_to_here();

					for (int i = 0; i != temp_jump_vector_to_start.size(); i++) {
						//追加するjumpの始点にちゃんとプロトンがいれば追加
						if (vector_finder(Diffusionspecie::diffusion_siteid_now_list, temp_jump_vector_to_start[i].get_start_site_id())) {
							jumps_possible.push_back(temp_jump_vector_to_start[i]);
						}
					}

					//end_wasからのjumpsを追加
					vector<Jump> temp_jump_vector_from_end = sites[end_was-1].get_jumps_from_here();

					for (int i = 0; i != temp_jump_vector_from_end.size(); i++) {
						//追加するjumpの終点にプロトンがいない、かつ終点がstart_wasじゃなければ追加
						if (!vector_finder(Diffusionspecie::diffusion_siteid_now_list, temp_jump_vector_from_end[i].get_end_site_id()) 
							&&
							temp_jump_vector_from_end[i].get_end_site_id() != start_was) {
							jumps_possible.push_back(temp_jump_vector_from_end[i]);
						}
					}
			

				/*	//jumps_impossible_tmpをjumps_possibleに追加する
					for (int i = 0; i != jumps_impossible_tmp.size(); i++) {
						jumps_possible.push_back(jumps_impossible_tmp[i]);
					}

					//jumps_impossible_tmpを空にする
					jumps_impossible_tmp.clear();
					jumps_impossible_tmp.shrink_to_fit();
				*/
			
				}



				
				
			}


			//起きうるジャンプについて、freqの総和をとりfreq_sumに格納する
			
			double freq_sum = 0.0;
			for (int i = 0, n = jumps_possible.size(); i != n; i++) {
				freq_sum += jumps_possible[i].get_freq();
			}

          //確認用
/*
			cout << "\t" << "diffusion_siteid_now_list = " ;
			for (int i = 0; i != Diffusionspecie::diffusion_siteid_now_list.size(); i++) {
				cout << "\t" << Diffusionspecie::diffusion_siteid_now_list[i] ; 
			}
			cout << endl;
			//確認用
			cout << "blocking_list = " ;
			for (auto itr = Diffusionspecie::blocking_list.begin(); itr != Diffusionspecie::blocking_list.end() ; itr++) {
				cout << *itr << "," ;
			}
			cout << endl;
			for (int i = 0; i != jumps_possible.size(); i++) {
				cout << "\t" << "\t" << "start = " << jumps_possible[i].get_start_site_id() ; 
				cout << "\t" << "\t" << "end = " << jumps_possible[i].get_end_site_id() ; 
				cout << "\t" << "\t" << "jumps_possible[" << i << "] = " << jumps_possible[i].get_freq() << endl;
			}
			cout << endl;

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
			for (int l = 0, n = jumps_possible.size() ; l != n; l++) {
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
			//diffusion_species[x].sum_squared_distanceにAng.に直した変位の2乗を追加
			Eigen::Vector3d eigen_happend_jump_vector;
			eigen_happend_jump_vector = transcoords(jumps_possible[jump_happen_number].get_jump_vector(),lattice_matrix);
			vector<double> cartesian_happend_jump_vector = eigen2vector(eigen_happend_jump_vector);

			for (int i = 0, n = diffusion_species.size(); i != n ; i++) {
				if ( diffusion_species[i].get_diffusion_id() == diff_tmp ) {
					vector<double> jump_total_tmp(3,0.0);
					vector<double> sum_squared_distance_tmp(3,0.0);
					for (int k = 0; k != jump_total_tmp.size(); k++) {
						jump_total_tmp[k] = diffusion_species[i].get_jump_total()[k] + jumps_possible[jump_happen_number].get_jump_vector()[k];
						sum_squared_distance_tmp[k] = diffusion_species[i].get_sum_squared_distance()[k] + pow(cartesian_happend_jump_vector[k],2);
						
					}
					diffusion_species[i].set_jump_total(jump_total_tmp);
					diffusion_species[i].set_sum_squared_distance(sum_squared_distance_tmp);

					//diffusion_counterを1つ増加
					diffusion_species[i].diffusion_counter++;


					//cout << "diffusion_species[" << i << "].diffusion_counter = " << diffusion_species[i].diffusion_counter << endl;
					


				}

			}




/*			for (int i = 0; i != Diffusionspecie::diffusion_siteid_now_list.size() ; i++) {

				cout << "diffusion_siteid_now_list[" << i  << "] = " << Diffusionspecie::diffusion_siteid_now_list[i] << endl;
			}
*/


			//diffusion_siteid_now_listを更新する
			vector<int>::iterator itr;
			start_was = jumps_start_id;
			end_was = jumps_end_id;
			itr = find(Diffusionspecie::diffusion_siteid_now_list.begin(), Diffusionspecie::diffusion_siteid_now_list.end(), start_was);
			if (itr == Diffusionspecie::diffusion_siteid_now_list.end()) cout << "search failed" << endl;
			int wanted_index = distance(Diffusionspecie::diffusion_siteid_now_list.begin(), itr);
			//始点と終点のidを交換することで更新完了
			Diffusionspecie::diffusion_siteid_now_list[wanted_index] = jumps_end_id;

			//blocking_listを更新する
			if (blocking_yes) {

				//start_wasのmate_listを削除する
				//cout << '\t' << "start_was = " << start_was << endl;
				//cout << '\t' << "start_was_mate_list = " ;
				for (int l = 0, n = sites[start_was-1].get_blocking_mate_list().size(); l != n; l++) {
					//cout << sites[start_was-1].get_blocking_mate_list()[l] << "," ;
					Diffusionspecie::blocking_list.erase(sites[start_was-1].get_blocking_mate_list()[l]);
				}
				//cout << endl;

				//end_wasのmate_listを追加する
				//cout << '\t' << "end_was = " << end_was << endl;
				//cout << '\t' <<  "end_was_mate_list = " ;
				for (int l = 0, n = sites[end_was-1].get_blocking_mate_list().size(); l != n; l++) {
					//cout << sites[end_was-1].get_blocking_mate_list()[l] << "," ;
					Diffusionspecie::blocking_list.insert(sites[end_was-1].get_blocking_mate_list()[l]);
				}
				//cout << endl;
			}
		
			


			//時間を更新する
			double rho_2 = random0to1(mt);
			//cout << "\t" << "\t" << "rho_2 = " << rho_2 << endl;

			double dt = -log(rho_2)/freq_sum;
			total_time += dt;
			//cout << "dt = " << dt << endl;
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

			//cout << "loop_counter " << loop_counter << " finished" << endl;

			




		}

		//ループ終了後
		//cout << "\t" << "loop_counter finished" << endl;
		//cout << endl;
		
		//cout << "\t" << "total_time = " << scientific << total_time << defaultfloat << endl;
		
		//拡散係数を出力およびファイルに出力する


		//自己拡散係数を求めるために拡散種ごとのjump_totalを準備する
		vector<double> jump_total_all(3,0.0);

		//拡散種1つ1つに対し操作を行う
		for (int j = 0, n = diffusion_species.size(); j != n; j++) {


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

		//平均変位を計算する
		vector<double> average_displacement(3,0.0);

		for (int j = 0; j != jump_total_all.size(); j++) {
			average_displacement[j] = jump_total_all[j] / diffusion_species.size();
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

		//transcoords関数でaverage_displacementを分率座標からcartesian座標に直す
		Eigen::Vector3d average_displacement_eigen;
		average_displacement_eigen = transcoords(average_displacement, lattice_matrix);

		//全原子の平均変位を出力するために、std::vectorに変換して格納しておく
		vector<double> average_displacement_cartesian = eigen2vector(average_displacement_eigen);
		average_displacement_vector.push_back(average_displacement_cartesian);


		//電場ありの場合、伝導度テンソルを出力
		vector<double> Sigma_x(3,0.0);
		if (E_field_yes) {
			double concentration = diffusion_species.size() / lattice_matrix.determinant();
			//cout << "lattice_matrix.determinant() = " << lattice_matrix.determinant() << endl;
			//cout << "concentration = " << concentration << endl;
		
			for (int j = 0; j != Sigma_x.size(); j++) {
				
				Sigma_x[j] = average_displacement_eigen[j];
				Sigma_x[j] *= q_charge * concentration / (E_field_strength * total_time);
				Sigma_x[j] *= pow(10,8); //Å^-1をcm^-1に変換
				//cout << "Sigma_x[" << j << "]  [S/cm] = " << Sigma_x[j] << endl;


			}
			
			Sigma_vector.push_back(Sigma_x);


		}



		//電場なしの場合、自己拡散係数を出力
		else {
			displacement_vector *= pow(10,-8); //[Ang.]を[cm]に変換

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



		//それぞれの拡散種について、変位を出力する
		for (int i = 0, n = diffusion_species.size(); i != n ; i++) {
			
			//the_number_of_KMCを出力する
			ofs_ave_dis << step_counter << ",";

			//diffusion_idを出力する
			ofs_ave_dis << diffusion_species[i].get_diffusion_id() << "," ;
			
			//transcoords関数でjump_totalを分率座標からcartesian座標に直す
			Eigen::Vector3d displacement_vector;
			displacement_vector = transcoords(diffusion_species[i].get_jump_total(),lattice_matrix);

			for (int j = 0; j != diffusion_species[i].get_jump_total().size(); j++) {
				ofs_ave_dis << fixed << setprecision(10) << displacement_vector(j) << defaultfloat << "," ;
				
			}

			for (int j = 0; j != diffusion_species[i].get_sum_squared_distance().size(); j++) {
				ofs_ave_dis << fixed << setprecision(10) << diffusion_species[i].get_sum_squared_distance()[j] << defaultfloat << "," ;
				
			}

			

			//diffusion_counterを出力する
			ofs_ave_dis << diffusion_species[i].diffusion_counter ;

			//１つの拡散種の変位を出力し終わったら改行する
			ofs_ave_dis << endl;

		}
		ofs_ave_dis << "##############################"<< endl;

		//最初だけ代入する定数値
		if (step_counter == 1) {
			map_for_output["concentration"] = diffusion_species.size()/lattice_matrix.determinant();
			map_for_output["temperture"] = temperture;
			map_for_output["ion_charge"] = ion_charge;

			if (E_field_yes) { 
				map_for_output["Efield_strength_Ang"] = E_field_strength;
				map_for_output["Efield_strength_cm"] = E_field_strength*pow(10,8);
			}

		}

		//今回のstepのトータルシミュレーション時間を記録しておく
		total_time_for_output.push_back(total_time);
		

		//sitesをリセットする
		for (int i = 0; i != sites.size(); i++ ) {
			sites[i].set_site_atom(1);
			sites[i].set_diffusion_id(-1);
		}
	
		//diffusion_siteid_now_listをリセットする
		Diffusionspecie::diffusion_siteid_now_list.clear();

		//blocking_listをリセットする
		Diffusionspecie::blocking_list.clear();



		//実行時間の計測
	//	end = chrono::system_clock::now();

	//	auto time = end - start;

	//	time_stamp = chrono::system_clock::to_time_t(start);
		//cout <<  ctime(&time_stamp);

	//	auto msec = chrono::duration_cast<chrono::microseconds>(time).count();
		//cout <<  msec << " msec" << endl;
		//cout << endl;

		//KMCが何回終わったか
		//cout << '\t' <<  step_counter << " times KMC finished" << endl;


	}
	
	//OUTPUTファイルをまとめる
	ofs_output << "#This is OUTPUT file written by toml format." << endl;
	//ofs_output << "#KMC = " << step_counter << " times" << endl;
	ofs_output << endl;
	ofs_output << "#total_time [s] (average of all KMCs)" << endl;
	ofs_output << "total_time = " << scientific << accumulate(total_time_for_output.begin(), total_time_for_output.end(), 0.0) / total_time_for_output.size()  << endl; //シミュレーション回数分の平均値を出す
	ofs_output << endl;
	ofs_output << "#concentration [/Ang.^3] " << endl;
	ofs_output << "concentration = " << map_for_output["concentration"] <<  endl;
	if (E_field_yes) {
		ofs_output << endl;
		ofs_output << "#Efield_strength [V/Ang.] " << endl;
		ofs_output << "Efield_strength.Ang = " << map_for_output["Efield_strength_Ang"] << endl;
		ofs_output << endl;
		ofs_output << "#Efield_strength [V/cm] "  << endl;
		ofs_output << "Efield_strength.cm = " << map_for_output["Efield_strength_cm"] << endl;
	}
	ofs_output << defaultfloat;
	ofs_output << endl;
	ofs_output << "#temperture [K] " << endl;
	ofs_output << "temperture = " << map_for_output["temperture"] << endl;
	ofs_output << endl;
	ofs_output << "#ion_charge 拡散種の電荷(整数) " << endl;
	ofs_output << "ion_charge = " << map_for_output["ion_charge"] << endl;
	ofs_output << endl;
	ofs_output << "#dimensionality : the dimension in this diffusion" << endl;
	ofs_output << "dimensionality = " << dimensionality << endl;
	ofs_output << endl;

	ofs_output << fixed << setprecision(10) ;



	//電場なしの場合
	if (!E_field_yes) {


		//トレーサー拡散係数
		ofs_diff << "tracer diffusion coefficient [cm^2/s]" << endl;
		double D_t_x = 0;
		double D_t_y = 0;
		double D_t_z = 0;
		for (int i = 0; i != D_t_3d_vector.size(); i++) {

			D_t_x += D_t_3d_vector[i][0];
			D_t_y += D_t_3d_vector[i][1];
			D_t_z += D_t_3d_vector[i][2];
				
			if (i+1 == D_t_3d_vector.size()) {

				D_t_x /= D_t_3d_vector.size();
				D_t_y /= D_t_3d_vector.size();
				D_t_z /= D_t_3d_vector.size();
			}

		}

		ofs_diff << "Dx = " << D_t_x << endl; 
		ofs_diff << "Dy = " << D_t_y << endl; 
		ofs_diff << "Dz = " << D_t_z << endl; 
		ofs_diff << endl;

		ofs_output << "#tracer diffusion coefficient [cm^2/s]" << endl;
		ofs_output << "tracer_D = " << "[" << D_t_x << "," << D_t_y << "," << D_t_z << "]" << endl;
		ofs_output << endl;
		
		
		//自己拡散係数
		ofs_diff << "self diffusion coefficient [cm^2/s]" << endl;
		double D_j_x = 0;
		double D_j_y = 0;
		double D_j_z = 0;
		for (int i = 0; i != D_j_3d_vector.size(); i++) {

			D_j_x += D_j_3d_vector[i][0];
			D_j_y += D_j_3d_vector[i][1];
			D_j_z += D_j_3d_vector[i][2];
				
			if (i+1 == D_j_3d_vector.size()) {

				D_j_x /= D_j_3d_vector.size();
				D_j_y /= D_j_3d_vector.size();
				D_j_z /= D_j_3d_vector.size();
			}

		}

		ofs_diff << "Dx = " << D_j_x << endl; 
		ofs_diff << "Dy = " << D_j_y << endl; 
		ofs_diff << "Dz = " << D_j_z << endl; 

		ofs_output << "#self diffusion coefficient [cm^2/s]" << endl;
		ofs_output << "self_D = " << "[" << D_j_x << "," << D_j_y << "," << D_j_z << "]" << endl;
		ofs_output << endl;


		
		//トレーサー拡散係数からトレーサー伝導度
		std::vector<double> tracer_Sigma(3, 0.0);
		tracer_Sigma[0] = NernstEinstein_DtoSigma(D_t_x, concentration, temperture, ion_charge, dimensionality);
		tracer_Sigma[1] = NernstEinstein_DtoSigma(D_t_y, concentration, temperture, ion_charge, dimensionality);
		tracer_Sigma[2] = NernstEinstein_DtoSigma(D_t_z, concentration, temperture, ion_charge, dimensionality);


		ofs_sigma << "tracer ionic conductivity " << conductivity_unit << endl;
		ofs_sigma << "Sigma_x = " << tracer_Sigma[0] << endl;
		ofs_sigma << "Sigma_y = " << tracer_Sigma[1] << endl;
		ofs_sigma << "Sigma_z = " << tracer_Sigma[2] << endl;
		ofs_sigma << endl;

		ofs_output << "#tracer ionic conductivity " << conductivity_unit << endl;
		ofs_output << "tracer_Sigma = " << "[" << tracer_Sigma[0] << "," << tracer_Sigma[1] << "," << tracer_Sigma[2] << "]" << endl;
		ofs_output << endl;

		//自己拡散係数から自己伝導度
		std::vector<double> self_Sigma (3,0.0);
		self_Sigma[0] = NernstEinstein_DtoSigma(D_j_x, concentration, temperture, ion_charge, dimensionality);
		self_Sigma[1] = NernstEinstein_DtoSigma(D_j_y, concentration, temperture, ion_charge, dimensionality);
		self_Sigma[2] = NernstEinstein_DtoSigma(D_j_z, concentration, temperture, ion_charge, dimensionality);
		ofs_sigma << "self ionic conductivity " << conductivity_unit << endl;
		ofs_sigma << "Sigma_x = " << self_Sigma[0] << endl;
		ofs_sigma << "Sigma_y = " << self_Sigma[1] << endl;
		ofs_sigma << "Sigma_z = " << self_Sigma[2] << endl;
		ofs_sigma << endl;

		ofs_output << "#self ionic conductivity " << conductivity_unit << endl;
		ofs_output << "self_Sigma = " << "[" << self_Sigma[0] << "," << self_Sigma[1] << "," << self_Sigma[2] << "]" << endl;
		ofs_output << endl;
	}



	//電場ありの場合
	if (E_field_yes) {


		//ofs_diff << scientific << "concentration = " << concentration << endl;
		//ofs_diff << scientific << "temperture = " << temperture << endl;
		//ofs_diff << scientific << "ion_charge = " << ion_charge << endl;


/*		//電場の向きにより出力する伝導度を変える(main関数の最初に書いた)
		string axis;
		switch (E_field_axis) {
			case 1 : axis = "x"; break;
			case 2 : axis = "y"; break;
			case 3 : axis = "z"; break;
			case -1 : axis = "-x"; break;
			case -2 : axis = "-y"; break;
			case -3 : axis = "-z"; break;
		}
*/

		//伝導度x成分
		ofs_sigma << "chemical ionic conductivity " << conductivity_unit << endl;
		ofs_diff << "chemical diffusion coefficient [cm^2/s]" << endl;

		vector<double> Sigma_vector_total(3,0.0);
		string j_to_axis;
		for (int j = 0; j != Sigma_vector_total.size(); j++) {
			switch (j) {
				case 0 : j_to_axis = "x"; break;
				case 1 : j_to_axis = "y"; break;
				case 2 : j_to_axis = "z"; break;
			}
			for (int i = 0; i != Sigma_vector.size(); i++) {

				Sigma_vector_total[j] += Sigma_vector[i][j];

				if (i+1 == Sigma_vector.size()) {
					Sigma_vector_total[j] /= Sigma_vector.size();
					ofs_sigma << "Sigma_" << axis <<  j_to_axis << " = " << scientific << Sigma_vector_total[j] << endl;
					ofs_diff << "D" << axis <<  j_to_axis << " = " << scientific << NernstEinstein_SigmatoD(Sigma_vector_total[j], concentration, temperture, ion_charge, dimensionality)  << endl;

				}
			
			}

			/*double variance_sigma = 0.0 ;
			double standard_deviation = 0.0 ;
			for (int i = 0; i != Sigma_vector.size(); i++) {
				variance_sigma += pow(Sigma_vector[i][j] - Sigma_vector_total[j], 2);

				if (i+1 == Sigma_vector.size()) {
					variance_sigma /= Sigma_vector.size();
					standard_deviation = sqrt(variance_sigma);
					ofs_sigma << "Standard deviation = " << scientific << standard_deviation << endl;
					
				}

			}
			*/
		}

		ofs_output << "#chemical ionic conductivity " << conductivity_unit << endl;
		ofs_output << "chemical_Sigma = " << "[" << Sigma_vector_total[0] << "," << Sigma_vector_total[1] << "," << Sigma_vector_total[2] << "]" << endl;
		ofs_output << endl;
		ofs_output << "#chemical diffusion coefficient [cm^2/s]" << endl;
		ofs_output << "chemical_D = " << "[" << NernstEinstein_SigmatoD(Sigma_vector_total[0], concentration, temperture, ion_charge, dimensionality) << 
 "," << NernstEinstein_SigmatoD(Sigma_vector_total[1], concentration, temperture, ion_charge, dimensionality) << "," << NernstEinstein_SigmatoD(Sigma_vector_total[2], concentration, temperture, ion_charge, dimensionality) << "]" << endl;
		ofs_output << endl;
		//化学拡散係数
		//ofs_diff << "chemical diffusion coefficient (cm^2/s)" << endl;
		

		
		



	}

	/*if (E_field_yes) {
		ofs_ave_dis << '\t' << "E_field_strength = " << E_field_strength << endl;
		ofs_ave_dis << '\t' << "E_field_axis = " <<  E_field_axis << " (+x=1, +y=2, +z=3, -x=-1, -y=-2, -z=-3)" << endl;
	}
	else {
		ofs_ave_dis << '\t' << "No E_field" << endl;
	}

*/

	//全ての原子の平均変位を出力しておく
	vector<double> av_dis_sum(3,0.0);
	for (int i = 0; i != av_dis_sum.size(); i++) {
		for (int j = 0; j != average_displacement_vector.size(); j++) {
			av_dis_sum[i] += average_displacement_vector[j][i];
		}
		av_dis_sum[i] /= average_displacement_vector.size();
		
	}
	

	ofs_output << "#mean displacement [Ang.] " << endl;
	ofs_output << "mean_displacement = " << "[" << av_dis_sum[0] << "," << av_dis_sum[1] << "," << av_dis_sum[2] << "]" << endl;
	ofs_output << endl;
	ofs_output << defaultfloat ;


	
	//実行時間の計測
	
	//プログラム開始時刻を表示
	end = chrono::system_clock::now();
	
	//開始時間のmsecを求める
	end_msec = chrono::duration_cast<chrono::milliseconds>(end.time_since_epoch());
	all_msec = end_msec.count();
	msec = all_msec % 1000 ;

	time_stamp = chrono::system_clock::to_time_t(end);
    timer = localtime(&time_stamp);
	cout << endl;
	cout << "#################################" << endl;
	cout << timer->tm_year + 1900 << "-" 
		<< timer->tm_mon +1 << "-" 
		<< timer->tm_mday << " " 
		<< timer->tm_hour << ":"
		<< timer->tm_min << ":"
		<< setfill('0') << right << setw(2) << timer->tm_sec << setfill(' ') << "." 
		<< setfill('0') << right << setw(3) << msec << setfill(' ') << endl;
	cout << "Program finished!" << endl;
	auto time = end - start;
	auto sec = chrono::duration_cast<chrono::seconds>(time).count();
	cout << "Execution time = " << sec << " sec" << endl;
	cout << "#################################" << endl;
	cout << endl;
}

