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
	int diffusion_start_site;
	int diffusion_end_site;
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

//Eigen::VectorXdをstd::vectorに型変換する関数
vector<double> eigen2vector(Eigen::Vector3d eigen_vector) 
{
	vector<double> vector_cartesian(3,0.0);
	for (int i = 0 ; i != vector_cartesian.size() ; i++) {
		vector_cartesian[i] = eigen_vector(i) ; 
	}

	return vector_cartesian;
}

//ある値がvector内の要素に含まれているか否か判定する関数
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
	//tomlファイルを読み込めるか確認しておく
	//auto toml_file = toml::parse("INPUT");
	//int MCSP = toml::find<int>(toml_file,"MCSP");

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
	ofs_ave_dis << "#the number of KMC, diffusion_id, displacement in x direction [Ang.], displacement in y direction [Ang.], displacement in z direction [Ang.], sum of squared displacement of each jumps in x, y, z [Ang.^2], start_site, end_site, jump_counter of total, rotation, hopping [times]" << endl;
	ofs_ave_dis << "KMC_times,diffusion_id,dx,dy,dz,sum_x2,sum_y2,sum_z2,start_site,end_site,jump_counter,rot_counter,hop_counter" << endl;

	//結果を出力するアウトプットファイルを作成する
	ofstream ofs_diff("DiffusionCoefficient", ios::app);

	//結果を出力するアウトプットファイルを作成する
	ofstream ofs_sigma("IonicConductivity", ios::app);


	//make DiffusionCoefficient vector
	vector< vector<double> > D_t_3d_vector;
	vector< vector<double> > D_j_3d_vector;
	vector< vector<double> > D_c_3d_vector;

	//make ElectricalConductivity vector
	vector< vector<double> > Sigma_vector;

	//make average_displacement vector
	vector< vector<double> > average_displacement_vector;

	//tomlファイルとしてINPUTを読み込む
	auto toml_file = toml::parse("INPUT");
	//INPUTを読み込む, find_or関数は該当するパラメータ(第2引数)が存在しない場合, デフォルト値(第3引数)を読み込む
	long long mcsp = toml::find_or<int>(toml_file,"MCSP", 0);
	int average = toml::find_or<int>(toml_file,"AVERAGE", 0);
	long long p_place_n = toml::find_or<int>(toml_file,"NDIFFS", 0);
	int E_field_yes = toml::find_or<int>(toml_file,"EFIELDON", 0);
	double correct_constant_for_pow = toml::find_or<double>(toml_file,"CORRECT", 0);
	int temperture = toml::find_or<int>(toml_file,"TEMP", 0);
	int E_field_axis = toml::find_or<int>(toml_file,"AXIS", 0);
	double distance_jump = toml::find_or<double>(toml_file,"DISTANCEJUMP", 1); //単位は[Å]
	int site_PE_read_yes = toml::find_or<int>(toml_file,"SITEPEREAD", 0);
	int blocking_list_read_yes = toml::find_or<int>(toml_file,"BLOCKINGLISTREAD", 0);
	int blocking_yes = toml::find_or<int>(toml_file,"BLOCKING", 0);
	int rot_hop_count_yes = toml::find_or<int>(toml_file,"ROTHOPCOUNT", 0);
	int dimensionality = toml::find_or<int>(toml_file,"DIMENSIONALITY", 3);
	auto anti_drift_Efield = toml::find<std::vector<double>>(toml_file,"ANTIDRIFT");

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
	//cout << "typeid(toml_file) = " << typeid(toml_file).name() << endl;
	if (!typeid(toml_file).name()) {
		cerr << "Could not find file INPUT" << endl;
		abort();
	}

	//return 0;

	cout << "INPUT read" << endl;
	cout << "\t" << "Monte Carlo Step per Particle : MCSP = " << mcsp << endl;
	cout << "\t" << "the number of particles for calculate ensemble average :  AVERAGE = " << average << endl;
	cout << "\t" << "the number of executing KMCs : N_KMC = " << step_max << endl;
	cout << "\t" << "Monte Carlo Step (=MCSP*NDIFFS) :  MCS = " << loop_max << endl;
	cout << "\t" << "the number of placing diffusion atoms :  NDIFFS = " << p_place_n << endl;
	cout << "\t" << "Electrical field strength :  EFIELD = " << scientific << E_field_strength << " [V/Å] " << endl;
	cout << "\t" << "Electrical field strength :  EFIELD = " << E_field_strength*pow(10,8) << defaultfloat << " [V/cm] " << endl;
	cout << "\t" << "temperture : TEMP = " << temperture << endl;
	cout << "\t" << "the correct_constant for adjusting jump frequency gamma (default = 0.1) :  correct_constant = " << correct_constant << endl;
	cout << "\t" << "the direction of E_field : E_field_axis = " << E_field_axis << " (+x=1, +y=2, +z=3, -x=-1, -y=-2, -z=-3)" << endl;
	cout << "\t" << "whether sitePE is read or not = " << site_PE_read_yes << " (1=read、0=not read)"<< endl;
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
			stringstream ss_line; 
			ss_line << line;
			string s;
			int s_counter = 0;

			//座標をコピーするためのfrac_dblvecを作成する
			vector<double> frac_dblvec(3,0.0);
			for (int i = 0; i <= 2; i++) {
				ss_line >> s;
				frac_dblvec[i] = stod(s);
			}

			//siteのidは「現在の行数-"DIRECT"の行数」
			int site_id_tmp = n_lines - DIRECT_num;

			//sites自体は0から始まるので、1つずらして代入する
			sites[site_id_tmp-1].set_site_id(site_id_tmp);
			sites[site_id_tmp-1].set_site_frac_coords(frac_dblvec);
		}

		n_lines += 1;
		
	}
	
		//格子ベクトル確認用
		cout << "lattice_matrix = " << endl;
		for (int i = 0; i <= 2; i++) {
			cout << "\t"  << lattice_matrix.row(i) << endl;
		} 
		cout << endl;


	//sitePE.datを読み込む
	//まずはPOSCARのサイト数と一致しているかを確認
	int site_pe_dat_total_number = 0;
	if (site_PE_read_yes) {
		ifstream site_pe_dat("sitePE.dat");

		if (!site_pe_dat) {
			cerr << "Could not find file sitePE.dat" << endl;
			abort();
		}

		string for_line_reader_sitepe;
		while (getline(site_pe_dat, for_line_reader_sitepe)) {
			
			//コメント行はスルー
			if (for_line_reader_sitepe[0] == '#' || for_line_reader_sitepe[0] == '/') {
				continue;
			}

			else {
				site_pe_dat_total_number++;
			}

		}

		if (site_pe_dat_total_number == site_total_number) {
			cout << "sitePE.dat read" << endl;
			cout << endl;
		}

		else {
			cerr << "sitePE.dat does not MATCH POSCAR." << endl;
			abort();
		}

	}


	
	//2次元vectorに、サイト番号とエネルギー,およびボルツマン因子を格納する(初期配置のため)
	vector< vector<double> > sitePE_dat(site_pe_dat_total_number, vector<double>(3));
	if (site_PE_read_yes) {

		ifstream ifs_sitepe("sitePE.dat");
		int n_lines = 0;
		while (getline(ifs_sitepe,line)) {

			//コメント行はスルー
			if (line[0] == '#' || line[0] == '/') {
				continue;
			}

			istringstream is(line);
			double site_id ;
			double site_pe;
			double boltzmann_factor;
			is >> site_id >> site_pe;
			boltzmann_factor = exp(-site_pe/(kb/ec*temperture));

			sitePE_dat[n_lines][0] = site_id;
			sitePE_dat[n_lines][1] = site_pe;
			sitePE_dat[n_lines][2] = boltzmann_factor;

			n_lines++;
			
			
		}
	}


	//blocking_list.csvを読み込む
	int csv_total_number = 0;
	if (blocking_list_read_yes) {
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

	if (blocking_list_read_yes) {
		ifstream ifs3("blocking_list.csv");
		n_lines = 0;
		while (getline(ifs3, line)){

			vector<string> strvec = split(line, ',');

			for (int i = 0; i != strvec.size(); i++){
				blocking_list_csv[n_lines].push_back(stoi(strvec[i]));
			}
			
			n_lines += 1;
		}	
	
		//blocking_list_csvをもとに、各Siteのblocking_mate_listに追加していく
		for (auto itr = blocking_list_csv.begin(); itr != blocking_list_csv.end(); itr++) {
			for (auto itr2 = (*(itr)).begin(); itr2 != (*(itr)).end(); itr2++){
				sites[*itr2-1].set_blocking_mate_list(*itr);
			}
		}

	}
	
	//拡散種の配置数がオーバーしていないか確認する
	//blockingするとき, blocking_list_csvの行数が有効サイト数の上限ゆえに, 上限を上回ってないか確認
	if (blocking_yes) {
		if (p_place_n > csv_total_number -1) {
			cerr << "NDIFF > effective site total number. Reduce NDIFF." << endl;
			abort();
		}
	}
	
	//blockingしないとき
	else {
		if (p_place_n > site_total_number -1) {
			cerr << "NDIFF > site total number. Reduce NDIFF." << endl;
			abort();
		}

	}



	//Site.site_frac_coordsをもとに、Jump.jump_vectorを生成する
	for (int i = 0; i != jumps.size(); i++) {

		//始点と終点のsite_idを取得
		int start_site_id_tmp = jumps[i].get_start_site_id();
		int end_site_id_tmp = jumps[i].get_end_site_id();

		//始点と終点の分率座標を取得
		//sites[i]のsite_idはi+1なので、tmpから1を引いておく
		vector<double> start_frac_coords = sites[start_site_id_tmp-1].get_site_frac_coords();
		vector<double> end_frac_coords = sites[end_site_id_tmp-1].get_site_frac_coords();

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
			}
		}

		//jump_vector_tmpをjump_vectorに代入する
		jumps[i].set_jump_vector(jump_vector_tmp);

	}



	//電場がかかっていた場合、もしくは勾配を打ち消す場合、ジャンプ頻度を補正する
	if (E_field_yes || anti_drift_Efield[0]) {

		//cartesian座標軸方向の単位ベクトルを作成し、電場の大きさをかけて電場ベクトルとする
		//変数E_field_axisによって作成する単位ベクトルを変える(x=0, y=1, z=2)
		
		Eigen::Vector3d E_field_vector(3); 

		if (E_field_yes) {
			switch (E_field_axis) {
				case 1 : E_field_vector << 1,0,0 ; break;
				case 2 : E_field_vector << 0,1,0 ; break;
				case 3 : E_field_vector << 0,0,1 ; break;
				case -1 : E_field_vector << -1,0,0 ; break;
				case -2 : E_field_vector << 0,-1,0 ; break;
				case -3 : E_field_vector << 0,0,-1 ; break;
			}

			E_field_vector *= E_field_strength;
		}

		else if (anti_drift_Efield[0]) {
			for (int i = 0; i != 3; i++) {
				E_field_vector(i) = anti_drift_Efield[i];
			}
		}
		cout << "E_field_vector = [ " << E_field_vector(0) << ", " << E_field_vector(1) << ", " << E_field_vector(2)  << " ]" << endl;

		//生成したジャンプに対し操作を行っていく
		for (int i = 0; i != jumps.size(); i++) {
			
			//transcoords関数でジャンプベクトルをfracからcartesianに直す
			Eigen::Vector3d jump_vector_cartesian = transcoords(jumps[i].get_jump_vector(), lattice_matrix);
			
			//ジャンプベクトル方向の電場の大きさを計算する
			double E_j_dot = E_field_vector.dot(jump_vector_cartesian);
			double E_along_jump_strength = E_j_dot / jump_vector_cartesian.norm();
		
			//ΔEmigを求める、1*で良いのはプロトンのみなことに注意
			double delta_E_mig = q_charge * E_along_jump_strength * jump_vector_cartesian.norm()/2;

			//ジャンプ頻度を補正する
			double fixed_jump_freq = jumps[i].get_freq() * exp(delta_E_mig/(kb*temperture));
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
	
	//blocking_list_csvとJMPDATAに整合性があるかを確認する(このブロックはBZY中のプロトン拡散にのみ適用可能)
	if (blocking_list_read_yes) {
		int match_counter = 0;
		//cout << "size = " << sites[0].get_jumps_from_here().size() << endl;
		for (int i = 0; i != sites[0].get_jumps_from_here().size() ; i++) {
			//cout << "for = " << i << " times" << endl;
			//blocking_mate_listにjumps_from_hereのget_end_site_idが2つ以上含まれていれば(回転経路に相当)整合性あり, BZYのみ
			if (vector_finder(sites[0].get_blocking_mate_list(), sites[0].get_jumps_from_here()[i].get_end_site_id()))
				match_counter++;
		}

		//cout << "match_counter = " << match_counter << endl;

		if (match_counter < 2) {
			cout << "blocking_list_csv does not match JMPDATA." << endl;	
			cerr << "blocking_list_csv does not match JMPDATA." << endl;	
			abort();
		}
	}

	//濃度を定義、セル内の拡散粒子数/セルの体積より単位は[1/Ang.^3]
	concentration = p_place_n / lattice_matrix.determinant();
	


	//ここからKMCシミュレーション開始, 指定回数のシミュレーションを繰り返す
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


		//乱数を準備しておく
		random_device rnd;
		mt19937 mt(rnd());

		
		//初期配置を生成する

		//まずはプロトンを配置する数を決定する(INPUTで指定されているはず)
		//もしINPUTでプロトンの配置数が設定されておらず、デフォルトの0が採用されていた場合
		if (p_place_n == 0) {

			cout << "NDIFFS not defined. decide NDIFFS." << endl;
			
			//1から(site_total_number-1)のうち、一様ランダムに1つを決定=プロトン配置数
			uniform_int_distribution<int> rnd_p_place_n(1,site_total_number-1);
			p_place_n = rnd_p_place_n(mt);
			cout << "プロトン配置数" << p_place_n << endl;
		

		}

		//sitePE_datのdeep_copyを取っておく(初期配置を決めるごとに変更するので)
		vector<vector<double>> sitePE_dat_deepcopy = sitePE_dat;

		//次に、どのサイトにプロトンを配置するかを決める
		vector<int> proton_place_number_vector;

		//sitePE.datに基づいて、ボルツマン因子に従ってプロトンを配置するサイトを選ぶ場合
		if (site_PE_read_yes) {
			
			//プロトンの初期配置数に達するまで繰り返す
			while (proton_place_number_vector.size() < p_place_n) { 

				//sitePE_dat[i][2]の和(ボルツマン因子の総和:分配関数)を計算する
				double part_func = 0;
				for (auto itr = sitePE_dat_deepcopy.begin(); itr != sitePE_dat_deepcopy.end(); itr++) {
					part_func += (*(itr))[2];
				}
				
				//0から1の乱数を生成する
				uniform_real_distribution<double> random0to1(0,1);
				double rho_0 = random0to1(mt);
				
				//sitePE_dat中のボルツマン因子を順に足し上げていき、和がrho_0*part_funcを超えたときの整数lを取得する
				double part_func_tmp = 0.0;
				int over_partial_number;
				for (int l = 0, n = sitePE_dat_deepcopy.size() ; l != n; l++) {
					part_func_tmp += sitePE_dat_deepcopy[l][2];

					if (part_func_tmp > rho_0 * part_func) {
						over_partial_number = l;
						break;
					}
				}

				//初期配置を格納したベクトルに、サイト番号を追加
				int added_proton_site_number = sitePE_dat_deepcopy[over_partial_number][0];
				proton_place_number_vector.push_back(added_proton_site_number);

				//初期配置として追加したサイトは分配関数から除く
				sitePE_dat_deepcopy.erase(sitePE_dat_deepcopy.begin() + over_partial_number);

				//area blockingを考慮する場合
				if (blocking_yes) {
					//先ほど初期配置としたサイトのarea blockingに関わるサイトもsitePE_datから削除する
					vector<int> temp_area_blocking_list = sites[added_proton_site_number-1].get_blocking_mate_list();
					//area blockingに関わるサイトを1つずつ探索する
					for (auto itr_blk = temp_area_blocking_list.begin() ; itr_blk != temp_area_blocking_list.end() ; itr_blk++){
						//cout << "area blocking site = " << *itr_blk << endl;
						auto itr_hit = find_if(sitePE_dat_deepcopy.begin(), sitePE_dat_deepcopy.end(), [=](const auto& row) { return (row[0] == *itr_blk); });
						//見つからなかった場合は何もしない
						if (itr_hit == sitePE_dat_deepcopy.end()) {
							//cout << "not found" << endl;
						}
						//見つかった場合は2つの位置の距離を計算後、削除する
						else {
							//cout << "found" << endl;
							auto dist = distance(sitePE_dat_deepcopy.begin(), itr_hit);
							sitePE_dat_deepcopy.erase(sitePE_dat_deepcopy.begin() + dist);
							
						}
						
					}
					
				}
			}
		}

		


		//ボルツマン因子は考慮せず、ランダムにプロトンを配置するサイトを選ぶ場合
		else {
			//blockingありの場合
			if (blocking_yes) {

				//乱数を生成し、blocking_list_csv(vector<vector<int>>)をシャッフルする
				shuffle( blocking_list_csv.begin(), blocking_list_csv.end(), mt );
				
				//先頭から抽出する
				for (int i = 0; i != p_place_n; i++) {
					shuffle ( blocking_list_csv[i].begin(), blocking_list_csv[i].end(), mt);
					proton_place_number_vector.push_back(blocking_list_csv[i][0]);
				}

			}

			//blockingなしの場合
			else {
				//(サイトの数-1)を要素にもつvectorを生成
				proton_place_number_vector.resize(site_total_number-1);
				for (int i = 0; i != proton_place_number_vector.size(); i++) {
					proton_place_number_vector[i] = i + 1;
				}

				//乱数を生成し、vectorをシャッフルする
				shuffle( proton_place_number_vector.begin(), proton_place_number_vector.end(), mt );
				

				//先頭からプロトンを配置する数分(p_place_n)だけ抜き出しソートする
				proton_place_number_vector.resize(p_place_n);
				sort( proton_place_number_vector.begin(), proton_place_number_vector.end() );

			}
		}
		
		

		//Diffusionspecieクラスのvectorをつくる(数はプロトン配置数=p_place_n)
		vector<Diffusionspecie> diffusion_species(p_place_n);
		//それぞれのidを[i]にたいしてi+1で設定する(1からp_place_nまでdiffusion_idとして通し番号をふる)
		for (int i = 0; i != diffusion_species.size(); i++) {

			int diffusion_start_site;
			diffusion_species[i].set_diffusion_id(i+1);
			diffusion_species[i].diffusion_start_site = proton_place_number_vector[i] ;
			diffusion_species[i].diffusion_counter = 0;
			diffusion_species[i].rotation_counter = 0;
			diffusion_species[i].hopping_counter = 0;
		}
		

		//プロトンを配置するSiteクラスのsite_atomをnumber_proton(本プログラムの最初で2を割り当てている)に変更する
		//と同時に、diffusion_idを設定する
		
		int d_id = 1;

		for (int i = 0; i != proton_place_number_vector.size(); i++) {
			for (int k = 0; k != sites.size(); k++) {

				//sites[k]=Siteクラスのsite_idが、プロトンを配置するsite_idかどうかを判定
				if (sites[k].get_site_id() == proton_place_number_vector[i]) {
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

					break;
				}
				else {
					continue;
				}
			}
		}
		
		
		//シミュレーション時間total_timeを定義
		double total_time = 0.0;

		//ここから1回のKMCシミュレーション内で指定されたステップ数を繰り返していく
		cout << endl;
		
		int start_was = 0;
		int end_was = 0;
		vector<Jump> jumps_possible;

		string processing;
		for (long long loop_counter = 1; loop_counter <= loop_max; loop_counter++) { 			
			
			//切り上げしておよそ10%ごとに#を出力する, log_coutに追記される, pythonでいうtqdmを簡略に実装
			if (loop_counter % ((loop_max+10-1)/10) == 0) {
				processing += "#";
				cout  << "\t" << processing << endl;
				if (loop_counter == loop_max ) {
					cout << endl;
				}
			}

			//1ループ目のみ
			if (loop_counter == 1) {
						
				//系で起きうる事象(今回はジャンプ)を列挙し、jumps_possibleに入れていく

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

			}


			//2ループ目以降
			else {

				//jumps_possibleを更新する

				//blockingあり, これは一般的なarea blockingに適用可能か?(高橋はBZYにのみ適用確認)
				if (blocking_yes) {
					
					auto itr = jumps_possible.begin();
					while (itr != jumps_possible.end()) {

						//start_was(直前のジャンプの開始点だったサイト)からのjumpを削除
						if ((*itr).get_start_site_id() == start_was) {
							itr = jumps_possible.erase(itr);
						}

						//end_was(直前のジャンプの終着点だったサイト)およびそのblocking_mate_listに向かうjumpを削除

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

						//end_was(直前のジャンプの終着点だったサイト)からのjumpを追加する
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
						
						//end_wasからのjumpを追加する
						vector<Jump> temp_jump_vector_from_end = sites[end_was-1].get_jumps_from_here();
						for (auto itr_2 = temp_jump_vector_from_end.begin(); itr_2 != temp_jump_vector_from_end.end(); itr_2++) {
							
							//回転経路内なら追加
							if (vector_finder(sites[end_was-1].get_blocking_mate_list(), (*itr_2).get_end_site_id())) {
								jumps_possible.push_back(*itr_2);
							}

							//ホッピング経路の場合
							else {
								//終点がblocking_listにない、かつ終点がstart_wasじゃなけば追加
								if (!Diffusionspecie::blocking_list.count((*itr_2).get_end_site_id())
									&&
									(*itr_2).get_end_site_id() != start_was) {

									jumps_possible.push_back(*itr_2);
								}
							}
						}

						//start_wasを含むblocking areaへのjumpを追加する
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
			
				}



				
				
			}


			//起きうるジャンプについて、ジャンプ頻度freqの総和をとりfreq_sumに格納する
			
			double freq_sum = 0.0;
			for (int i = 0, n = jumps_possible.size(); i != n; i++) {
				freq_sum += jumps_possible[i].get_freq();
			}


			//実際に起こるイベントを決める
			//0から1の乱数を生成する
			uniform_real_distribution<double> random0to1(0,1);
			double rho_1 = random0to1(mt);
			
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

					//rot_hop_count is 1 (True)
					if (rot_hop_count_yes) {
										
						//rotation
						if (vector_finder(sites[jumps_start_id-1].get_blocking_mate_list(), jumps_end_id)) {
							diffusion_species[i].rotation_counter++;
						}
							
						//hopping
						else {
							diffusion_species[i].hopping_counter++;
						}

					}

				}

			}


			//diffusion_siteid_now_listを更新する
			vector<int>::iterator itr;
			start_was = jumps_start_id; //直前のジャンプの始点サイトid
			end_was = jumps_end_id;     //直前のジャンプの終点サイトid
			itr = find(Diffusionspecie::diffusion_siteid_now_list.begin(), Diffusionspecie::diffusion_siteid_now_list.end(), start_was);
			if (itr == Diffusionspecie::diffusion_siteid_now_list.end()) cout << "search failed" << endl;
			int wanted_index = distance(Diffusionspecie::diffusion_siteid_now_list.begin(), itr);
			//始点と終点のidを交換することで更新完了
			Diffusionspecie::diffusion_siteid_now_list[wanted_index] = jumps_end_id;

			//blocking_listを更新する
			if (blocking_yes) {

				//start_wasのmate_listを削除する
				for (int l = 0, n = sites[start_was-1].get_blocking_mate_list().size(); l != n; l++) {
					Diffusionspecie::blocking_list.erase(sites[start_was-1].get_blocking_mate_list()[l]);
				}

				//end_wasのmate_listを追加する
				for (int l = 0, n = sites[end_was-1].get_blocking_mate_list().size(); l != n; l++) {
					//cout << sites[end_was-1].get_blocking_mate_list()[l] << "," ;
					Diffusionspecie::blocking_list.insert(sites[end_was-1].get_blocking_mate_list()[l]);
				}

			}
		
			

			//時間を更新する
			double rho_2 = random0to1(mt);

			double dt = -log(rho_2)/freq_sum;
			total_time += dt;

			//ループの最後の時
			if (loop_counter == loop_max) { 
				//end_siteを設定する
				for (int i = 0; i != Diffusionspecie::diffusion_siteid_now_list.size() ; i++) {
					int tmp_diff_id = sites[Diffusionspecie::diffusion_siteid_now_list[i]-1].get_diffusion_id() ;
					diffusion_species[tmp_diff_id-1].diffusion_end_site = Diffusionspecie::diffusion_siteid_now_list[i] ;
				}
			}

		}

		//1度のKMCシミュレーションでの既定ステップ数が終了したあと
		
		//拡散係数を出力およびファイルに出力する(この出力された拡散係数はドリフトを考慮しておらず、平均変位によるセンタリングも行っていない点に注意)

		//集団拡散係数および伝導度拡散係数を求めるために拡散種ごとのjump_totalを準備する
		vector<double> jump_total_all(3,0.0);

		//拡散粒子1つ1つに対し操作を行う
		for (int j = 0, n = diffusion_species.size(); j != n; j++) {

			//集団拡散係数(および電場勾配化では伝導度)を求めるため
			//拡散種ごとのjump_totalを足して合計変位を出す
			for (int i = 0; i != jump_total_all.size(); i++) {
				jump_total_all[i] += diffusion_species[j].get_jump_total()[i];
			}

			if (!E_field_yes) {   //電場なしのとき 

			//トレーサー拡散係数を計算し出力
			vector<double> D_t_3d = diffusion_species[j].get_D(lattice_matrix,total_time);

			D_t_3d_vector.push_back(D_t_3d);

			}



		} 

		//平均変位を計算する
		vector<double> average_displacement(3,0.0);

		for (int j = 0; j != jump_total_all.size(); j++) {
			average_displacement[j] = jump_total_all[j] / diffusion_species.size();
		}

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

			for (int j = 0; j != Sigma_x.size(); j++) {
				
				Sigma_x[j] = average_displacement_eigen[j];
				Sigma_x[j] *= q_charge * concentration / (E_field_strength * total_time);
				Sigma_x[j] *= pow(10,8); //Å^-1をcm^-1に変換
			}
			
			Sigma_vector.push_back(Sigma_x);


		}



		//電場なしの場合、集団拡散係数を出力
		else {
			displacement_vector *= pow(10,-8); //[Ang.]を[cm]に変換

			//jump_total_allを2乗したあと拡散種の数で割って、集団拡散係数を出力する
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

			//start_site,end_siteを出力する
			ofs_ave_dis << diffusion_species[i].diffusion_start_site << "," ;
			ofs_ave_dis << diffusion_species[i].diffusion_end_site << "," ;


			//diffusion_counterを出力する
			ofs_ave_dis << diffusion_species[i].diffusion_counter << "," ;

			if (rot_hop_count_yes) {
				ofs_ave_dis << diffusion_species[i].rotation_counter << "," ;
				ofs_ave_dis << diffusion_species[i].hopping_counter ;
			}

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
		
		
		//集団拡散係数
		ofs_diff << "collective diffusion coefficient [cm^2/s]" << endl;
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

		ofs_output << "#collective diffusion coefficient [cm^2/s]" << endl;
		ofs_output << "collective_D = " << "[" << D_j_x << "," << D_j_y << "," << D_j_z << "]" << endl;
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

		//集団拡散係数から集団伝導度
		std::vector<double> collective_Sigma (3,0.0);
		collective_Sigma[0] = NernstEinstein_DtoSigma(D_j_x, concentration, temperture, ion_charge, dimensionality);
		collective_Sigma[1] = NernstEinstein_DtoSigma(D_j_y, concentration, temperture, ion_charge, dimensionality);
		collective_Sigma[2] = NernstEinstein_DtoSigma(D_j_z, concentration, temperture, ion_charge, dimensionality);
		ofs_sigma << "collective ionic conductivity " << conductivity_unit << endl;
		ofs_sigma << "Sigma_x = " << collective_Sigma[0] << endl;
		ofs_sigma << "Sigma_y = " << collective_Sigma[1] << endl;
		ofs_sigma << "Sigma_z = " << collective_Sigma[2] << endl;
		ofs_sigma << endl;

		ofs_output << "#collective ionic conductivity " << conductivity_unit << endl;
		ofs_output << "collective_Sigma = " << "[" << collective_Sigma[0] << "," << collective_Sigma[1] << "," << collective_Sigma[2] << "]" << endl;
		ofs_output << endl;
	}



	//電場ありの場合
	if (E_field_yes) {

		//伝導度の各軸成分
		ofs_sigma << "ionic conductivity from drift velocity" << conductivity_unit << endl;
		ofs_diff << "conductivity diffusion coefficient [cm^2/s]" << endl;

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

		ofs_output << "#ionic conductivity from drift velocity" << conductivity_unit << endl;
		ofs_output << "Sigma = " << "[" << Sigma_vector_total[0] << "," << Sigma_vector_total[1] << "," << Sigma_vector_total[2] << "]" << endl;
		ofs_output << endl;
		ofs_output << "#conductivity diffusion coefficient [cm^2/s]" << endl;
		ofs_output << "D_sigma = " << "[" << NernstEinstein_SigmatoD(Sigma_vector_total[0], concentration, temperture, ion_charge, dimensionality) << 
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
	
	//プログラム終了時刻を表示
	end = chrono::system_clock::now();
	
	//終了時間のmsecを求める
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

