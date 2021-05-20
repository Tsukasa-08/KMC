#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>
#include <boost/lexical_cast.hpp>
#include <eigen3/Core>
#include <eigen3/LU>

using namespace std;

// Function for reading a data file.
vector<string> split(string& input, char delimiter)
{
	istringstream stream(input);
	string field;
	vector<string> result;
	while (getline(stream, field, delimiter)) {
		result.push_back(field);
	};
	return result;
};


// Define node class.
class Node {
private:
	int step;
	int site_id;
	int root_id;
	double time;
	double prob;
	vector<double> jmpVec;
	map<int,double> occ;

public:
	Node(int st, int sid, int rid, double ti, double pr, vector<double> jv, map<int,double> oc) {
		   step=st, site_id=sid, root_id=rid,time=ti, prob=pr, jmpVec=jv, occ=oc ;
	};

	int get_step() { return step; };
	int get_site_id() { return site_id; };
	int get_root_id() { return root_id; };
	double get_time() { return time; };
	double get_prob() { return prob; };
	vector<double> get_jmpVec() { return jmpVec; };
	map<int,double> get_occ() { return occ; };

	/*
	void update_time(double dt) { time += dt; };
	void update_prob(double prob_crt) { prob *= prob_crt; };
	void update_jmpVec(vector<double> jmpVec_crt) {
		   jmpVec.at(0) += jmpVec_crt.at(0);
		   jmpVec.at(1) += jmpVec_crt.at(1);
		   jmpVec.at(2) += jmpVec_crt.at(2);
	}
	*/
};


// Function for estimating average time of ocsiratory jumps between two sites.
double meanTimeOsc(Eigen::Matrix2d T, Eigen::Matrix2d P) {
	Eigen::Matrix2d A = (Eigen::Matrix2d::Identity() - P).inverse();
	Eigen::RowVector2d p0;
	p0 << 1.0,0.0;
	Eigen::Vector2d q;
	q << 1.0-P(0,1), 1.0-P(1,0);
	double time = p0 * A * T * P * A * q;
	return time;
};


// Function for estimating average time of ocsiratory jumps between two sites.
Eigen::Vector2d pEscapeOsc(Eigen::Matrix2d T, Eigen::Matrix2d P) {
	Eigen::Matrix2d A = (Eigen::Matrix2d::Identity() - P).inverse();
	Eigen::RowVector2d p0;
	p0 << 1.0,0.0;
	Eigen::Vector2d q;
	q << 1.0-P(0,1), 1.0-P(1,0);
	Eigen::Vector2d pEscape = p0 * A;
	pEscape(0) = pEscape(0) * q(0);
	pEscape(1) = pEscape(1) * q(1);
	return pEscape;
};


// Function for Combining the same jumps.
map<string,double> mkjmpdata(list<Node> &nodes) {
	map<string,double> jmpdata;
	for (auto itr3 = nodes.begin(); itr3 != nodes.end(); ++itr3) {
		int i = (*itr3).get_root_id() + 1;
		int j = (*itr3).get_site_id() + 1;
		vector<double> jmpVec = (*itr3).get_jmpVec();
		double time = (*itr3).get_time();
		double prob = (*itr3).get_prob();
		if (abs(jmpVec.at(0)) < 0.01 and abs(jmpVec.at(1)) < 0.01 and abs(jmpVec.at(2)) < 0.01) {
			;
		} else {
			stringstream ss;
			ss << fixed << i <<"," << j << ","
				 << jmpVec.at(0) << "," << jmpVec.at(1) << "," << jmpVec.at(2) << ",";
			if (jmpdata.count(ss.str())) {
				jmpdata[ss.str()] += prob/time;
			} else {
				jmpdata[ss.str()] = prob/time;
			};
		};
	};
	return jmpdata;
};


// Output files.
void output(map<string,double> &jmpData, string name) {
	ofstream outputfile(name);
	for (auto itr4 = jmpData.begin(); itr4 != jmpData.end(); ++itr4) {
		string key = itr4->first;
		double val = itr4->second;
		outputfile << key << scientific << setprecision(10) << val << "\n";
	};
	outputfile.close();
};


int main()
{
	// Read INPUT.
	int n_sites;
	int n_jmps;
	int n_loops;
	int sid_i;
	string name_jmpdata;
	string name_occ;

	ifstream ifs1("INPUT");
	string line;
	while (getline(ifs1, line)) {
		vector<string> strvec = split(line, ',');
		if ( strvec.at(0).find('#') != string::npos ) continue;
		n_sites = atoi(strvec.at(0).c_str());
		n_jmps = atoi(strvec.at(1).c_str());
		n_loops = atoi(strvec.at(2).c_str());
		sid_i = atoi(strvec.at(3).c_str());
		name_jmpdata = strvec.at(4);
		name_occ = strvec.at(5);
	}
	ifs1.close();

	
	// Read jumpdata.csv.
	vector< vector<int> > jmpSites_all(n_jmps,vector<int>(2,-1));
	vector< vector<double> > jmpVecs_all(n_jmps, vector<double>(3,0.0));
	vector<double> jmpFreqs_all(n_jmps,0.0);
	
	ifstream ifs2(name_jmpdata.c_str());
	int n_lines = 0;
	while (getline(ifs2, line)) {
		vector<string> strvec = split(line, ',');
		if ( strvec.at(0).find('#') != string::npos ) continue;
		jmpSites_all.at(n_lines).at(0) = atoi(strvec.at(0).c_str())-1;
		jmpSites_all.at(n_lines).at(1) = atoi(strvec.at(1).c_str())-1;
		jmpVecs_all.at(n_lines).at(0) = atof(strvec.at(2).c_str());
		jmpVecs_all.at(n_lines).at(1) = atof(strvec.at(3).c_str());
		jmpVecs_all.at(n_lines).at(2) = atof(strvec.at(4).c_str());
		jmpFreqs_all.at(n_lines) = atof(strvec.at(5).c_str());
		n_lines += 1;
	};
	ifs2.close();

	// Make jmpFreq matrix.
	vector< vector<double> > jmpFreqMat(n_sites, vector<double>(n_sites,0.0));
	for (int i = 0 ; i < jmpSites_all.size() ; i++) {
		jmpFreqMat.at(jmpSites_all.at(i).at(0)).at(jmpSites_all.at(i).at(1)) = jmpFreqs_all.at(i);
	};


	// Make jmpVec matrices (x,y,z).
	vector< vector<double> > jmpVecMat_x(n_sites, vector<double>(n_sites,0.0));
	vector< vector<double> > jmpVecMat_y(n_sites, vector<double>(n_sites,0.0));
	vector< vector<double> > jmpVecMat_z(n_sites, vector<double>(n_sites,0.0));
	for (int i = 0 ; i < jmpSites_all.size() ; i++) {
		jmpVecMat_x.at(jmpSites_all.at(i).at(0)).at(jmpSites_all.at(i).at(1)) = jmpVecs_all.at(i).at(0);
		jmpVecMat_y.at(jmpSites_all.at(i).at(0)).at(jmpSites_all.at(i).at(1)) = jmpVecs_all.at(i).at(1);
		jmpVecMat_z.at(jmpSites_all.at(i).at(0)).at(jmpSites_all.at(i).at(1)) = jmpVecs_all.at(i).at(2);
	};


	// Make adjacent vectors.
	vector< vector<int> > adj(n_sites);
	for (int i = 0 ; i < n_sites ; i++) {
		for (int j = 0 ; j < n_sites ; j++) {
			if ( jmpFreqMat.at(i).at(j) > 0.00001 ) {
				(adj.at(i)).push_back(j);
			};
		};
	};
	

	// Read occ_eq.csv.
	vector< vector<double> > occ_eq(n_sites,vector<double>(n_sites,0.0));
	
	ifstream ifs3(name_occ.c_str());
	n_lines = 0;
	while (getline(ifs3, line)) {
		vector<string> strvec = split(line, ',');
		if ( strvec.at(0).find('#') != string::npos ) continue;
		for (int i = 0 ; i < n_sites ; i++) {
			occ_eq.at(n_lines).at(i) = atof(strvec.at(i).c_str());
		}
		n_lines += 1;
	};
	ifs3.close();


	// Set initial sites.
	vector<int> sids_ini;
	if (sid_i == -1) {
		for (int i = 0 ; i < n_sites ; i++) {
			sids_ini.push_back(i);
		};
	} else {
		sids_ini.push_back(sid_i-1);
	};


	// Start main loops.
	// Initialize
	list<Node> nodes_prev = {};
	list<Node> nodes_crt = {};
	for (int i_ = 0 ; i_ < sids_ini.size() ; i_++) {
		int i = sids_ini.at(i_) ;
		map<int, double> occ_tmp;
		nodes_prev.push_back(Node(0, i, i, 0.0, 1.0, {0.0, 0.0, 0.0}, occ_tmp));
	};

	// 1st iteration (n = 1)
	int n_nodes = nodes_prev.size();
	cout << "1  " << n_nodes << "\n";
	for (int n_itr = 0 ; n_itr < n_nodes ; n_itr++) {
		cout << "n_itr: " << n_itr << "\n";
		int i = (nodes_prev.front()).get_site_id();
		int rid = (nodes_prev.front()).get_root_id();
		double t_prev = (nodes_prev.front()).get_time();
		double p_prev = (nodes_prev.front()).get_prob();
		vector<double> jmpVec_prev = (nodes_prev.front()).get_jmpVec();
		map<int,double> occ_prev = (nodes_prev.front()).get_occ();
		
		// Update dt and freqs
		double freq_sum = 0.0;
		vector<double> freqs((adj.at(i)).size(),0.0);
		for (int j_ = 0 ; j_ < (adj.at(i)).size() ; j_++ ) {
			int j = adj.at(i).at(j_);
			cout << "1st j: " << j << "  ";
			freqs.at(j_) = jmpFreqMat.at(i).at(j) * (1.0-occ_eq.at(i).at(j));
			cout << "hoge\n";
			freq_sum += freqs.at(j_);
		};
		cout << "\n";
		double dt = 1.0/freq_sum;
		map<int,double> occ_crt;
		occ_crt[i] = 0.0;
		for ( int j_ = 0 ; j_ < (adj.at(i)).size() ; j_++ ) {
			int j = adj.at(i).at(j_);
			cout << "2nd j: " << j << "  ";
			vector<double> jmpVec_crt = {jmpVec_prev.at(0)+jmpVecMat_x.at(i).at(j),
				                           jmpVec_prev.at(1)+jmpVecMat_y.at(i).at(j),
				                           jmpVec_prev.at(2)+jmpVecMat_z.at(i).at(j)};
			nodes_crt.push_back(Node(1,j,rid,t_prev+dt,p_prev*freqs[j_]*dt,jmpVec_crt,occ_crt));
		};
		cout << "\n";
		nodes_prev.pop_front();
	};
	nodes_prev = nodes_crt;
	nodes_crt = {};


	// Output jmpdata.csv.mod1
	map<string,double> jmpdata = mkjmpdata(nodes_prev);
	output(jmpdata, "jmpdata.csv.mod1");


	// 2nd Iteration (n = 2). Infinite oscilatory jumps between the two sites.
	n_nodes = nodes_prev.size();
	cout << "2  " << n_nodes << "\n";
	for (int n_itr = 0 ; n_itr < n_nodes ; n_itr++) {
		int i = (nodes_prev.front()).get_site_id();
		int rid = (nodes_prev.front()).get_root_id();
		double t_prev = (nodes_prev.front()).get_time();
		double p_prev = (nodes_prev.front()).get_prob();
		vector<double> jmpVec_prev = (nodes_prev.front()).get_jmpVec();

		// Determine dt and prob self-consistently.(1: occupied site i, 2 unoccupied site rid)
		vector<double> freqs1((adj.at(i)).size(),0.0);
		double dt1 = 0.0;
		double occ2 = 0.0;
		double p12 = 0.0;
		double ratio_max = 0.0;
		do {
			// Update occ
			double occ2_ave = 0.0;
			double inflow2 = 0.0;
			for (int k_ = 0 ; k_ < (adj.at(rid)).size() ; k_++) {
				int k = adj.at(rid).at(k_);
				if (k != i ) {
					inflow2 += occ_eq.at(i).at(k) * jmpFreqMat.at(k).at(rid);
				};
			};
			occ2 = occ_eq.at(i).at(rid)*(1.0-exp(-inflow2*dt1));
			occ2_ave = occ2/2.0;
			if (inflow2*dt1 > 1.0e-10) {
				occ2_ave = occ_eq.at(i).at(rid)*(1.0+(exp(-inflow2*dt1)-1.0)/(inflow2*dt1));
			};

			// Update dt and freqs
			double freq_sum = 0.0;
			ratio_max = 0.0;
			vector<double> freqs_prev = freqs1;
			for (int j_ = 0 ; j_ < (adj.at(i)).size() ; j_++ ) {
				int j = adj.at(i).at(j_);
				if (j == rid) {
					freqs1.at(j_) = jmpFreqMat.at(i).at(j) * (1.0-occ2_ave);
					p12 = freqs1.at(j_);
				} else {
					freqs1.at(j_) = jmpFreqMat.at(i).at(j) * (1.0-occ_eq.at(i).at(j));
				};
				freq_sum += freqs1.at(j_);
				double ratio = freqs1.at(j_)/freqs_prev.at(j_);
				if ( ratio < 1.0 ) {
					ratio = 1.0/ratio;
				};
				ratio_max = max(ratio_max,ratio);
			};
			dt1 = 1.0/freq_sum;
			p12 = p12*dt1;
		} while (ratio_max > 1.01);

		// Determine dt and prob self-consistently.(1: unoccupied site i, 2 occupied site rid)
		vector<double> freqs2((adj.at(rid)).size(),0.0);
		double dt2 = 0.0;
		double occ1 = 0.0;
		double p21 = 0.0;
		ratio_max = 0.0;
		do {
			// Update occ
			double occ1_ave = 0.0;
			double inflow1 = 0.0;
			for (int k_ = 0 ; k_ < (adj.at(i)).size() ; k_++) {
				int k = adj.at(i).at(k_);
				if (k != rid ) {
					inflow1 += occ_eq.at(rid).at(k) * jmpFreqMat.at(k).at(i);
				};
			};
			occ1 = occ_eq.at(rid).at(i)*(1.0-exp(-inflow1*dt2));
			occ1_ave = occ1/2.0;
			if (inflow1*dt2 > 1.0e-10) {
				occ1_ave = occ_eq.at(rid).at(i)*(1.0+(exp(-inflow1*dt2)-1.0)/(inflow1*dt2));
			};

			// Update dt and freqs
			double freq_sum = 0.0;
			ratio_max = 0.0;
			vector<double> freqs_prev = freqs2;
			for (int j_ = 0 ; j_ < (adj.at(rid)).size() ; j_++ ) {
				int j = adj.at(rid).at(j_);
				if (j == i) {
					freqs2.at(j_) = jmpFreqMat.at(rid).at(j) * (1.0-occ1_ave);
					p21 = freqs2.at(j_);
				} else {
					freqs2.at(j_) = jmpFreqMat.at(rid).at(j) * (1.0-occ_eq.at(rid).at(j));
				};
				freq_sum += freqs2.at(j_);
				double ratio = freqs2.at(j_)/freqs_prev.at(j_);
				if ( ratio < 1.0 ) {
					ratio = 1.0/ratio;
				};
				ratio_max = max(ratio_max,ratio);
			};
			dt2=1.0/freq_sum;
			p21 = p21*dt2;
		} while (ratio_max > 1.01);

		Eigen::Matrix2d T;
		T << dt1,0.0, 0.0,dt2;
		Eigen::Matrix2d P;
		P << 0.0,p12, p21,0.0;
		double meanT = meanTimeOsc(T,P);
		Eigen::Vector2d p_escape = pEscapeOsc(T,P);
		meanT = meanT + p_escape(0) * dt1 + p_escape(1) * dt2;

		// Add nodes adjacent to i
		for ( int j_ = 0 ; j_ < (adj.at(i)).size() ; j_++ ) {
			int j = adj.at(i).at(j_);
			float denom = (1.0-p12)/dt1;
			if (j != rid) {
				vector<double> jmpVec_crt = {jmpVec_prev.at(0)+jmpVecMat_x.at(i).at(j),
				                             jmpVec_prev.at(1)+jmpVecMat_y.at(i).at(j),
				                             jmpVec_prev.at(2)+jmpVecMat_z.at(i).at(j)};
				map<int,double> occ_crt;
				occ_crt[i] = 0.0;
				occ_crt[rid] = occ2;
				double freq_tmp = jmpFreqMat.at(i).at(j)*(1.0-occ_eq.at(i).at(j));
				nodes_crt.push_back(Node(2,j,rid,t_prev+meanT,p_prev*p_escape(0)*freq_tmp/denom,jmpVec_crt,occ_crt));
			};
		};

		// Add nodes adjacent to rid
		for ( int j_ = 0 ; j_ < (adj.at(rid)).size() ; j_++ ) {
			int j = adj.at(rid).at(j_);
			float denom = (1.0-p21)/dt2;
			if (j != i) {
				vector<double> jmpVec_crt = {jmpVecMat_x.at(rid).at(j),
				                             jmpVecMat_y.at(rid).at(j),
				                             jmpVecMat_z.at(rid).at(j)};
				map<int,double> occ_crt;
				occ_crt[rid] = 0.0;
				occ_crt[i] = occ1;
				double freq_tmp = jmpFreqMat.at(rid).at(j)*(1.0-occ_eq.at(rid).at(j));
				nodes_crt.push_back(Node(2,j,rid,t_prev+meanT,p_prev*p_escape(1)*freq_tmp/denom,jmpVec_crt,occ_crt));
			};
		};
		nodes_prev.pop_front();
	};
	nodes_prev = nodes_crt;
	nodes_crt = {};

	// Output jmpdata.csv.mod2
	jmpdata = mkjmpdata(nodes_prev);
	output(jmpdata, "jmpdata.csv.mod2");


	// 3rd Iteration and more (n: iteration, i: initial site, j: adjSites to i, k: adjSites to j)
	for (int n = 3; n <= n_loops ; n++) {
		int n_nodes = nodes_prev.size();
		cout << n << "  " << n_nodes  << "\n";
		for (int n_itr = 0 ; n_itr < n_nodes ; n_itr++) {
			int i = (nodes_prev.front()).get_site_id();
			int rid = (nodes_prev.front()).get_root_id();
			double t_prev = (nodes_prev.front()).get_time();
			double p_prev = (nodes_prev.front()).get_prob();
			vector<double> jmpVec_prev = (nodes_prev.front()).get_jmpVec();
			map<int,double> occ_prev = (nodes_prev.front()).get_occ();
			
			// Determine dt and prob self-consistently.
			vector<double> freqs((adj.at(i)).size(),0.0);
			map<int,double> occ_crt = occ_prev;
			double dt = 0.0;
			double ratio_max;
			do {
				// Update occ
				map<int,double> occ_ave;
				vector<int> del_occ;
				vector<int> del_occ_ave;
				for (auto itr2 = occ_prev.begin(); itr2 != occ_prev.end(); ++itr2) {
					int key = itr2->first;
					double val = itr2->second;
					if ( key != i and val < occ_eq.at(i).at(key)) {
						double inflow = 0.0;
						for (int k_ = 0 ; k_ < (adj.at(key)).size() ; k_++) {
							int k = adj.at(key).at(k_);
							if (k != i ) {
								if (occ_crt.count(k)) {
									inflow += occ_crt[k] * jmpFreqMat.at(k).at(key);
								} else {
									inflow += occ_eq.at(i).at(k) * jmpFreqMat.at(k).at(key);
								};
							};
							//double occ_noneq = val+inflow*dt;
							double occ_noneq = val+(occ_eq.at(i).at(key)-val)*(1.0-exp(-inflow*dt));
							double occ_ave_noneq = val+(occ_eq.at(i).at(key)-val)*(1.0-exp(-inflow*dt))/2.0;
							if (inflow*dt > 1.0e-10) {
								occ_ave_noneq = val+(occ_eq.at(i).at(key)-val)*(1.0+(exp(-inflow*dt)-1.0)/(inflow*dt));
							};
							if ( occ_noneq < occ_eq.at(i).at(key) ) {
								occ_crt[key] = occ_noneq;
							} else {
								del_occ.push_back(key);
							};
							if ( occ_ave_noneq < occ_eq.at(i).at(key) ) {
								occ_ave[key] = occ_ave_noneq;
							} else {
								del_occ_ave.push_back(key);
							};
						};
					} else {
						del_occ.push_back(key);
						del_occ_ave.push_back(key);
					};
				};
				for (int k = 0 ; k < del_occ.size() ; k++) {
					occ_crt.erase(del_occ.at(k));
				};
				for (int k = 0 ; k < del_occ_ave.size() ; k++) {
					occ_ave.erase(del_occ_ave.at(k));
				};

				// Update dt and freqs
				double freq_sum = 0.0;
				ratio_max = 0.0;
				vector<double> freqs_prev = freqs;
				for (int j_ = 0 ; j_ < (adj.at(i)).size() ; j_++ ) {
					int j = adj.at(i).at(j_);
					if (occ_ave.count(j) and occ_ave[j] < occ_eq.at(i).at(j)) {
						freqs.at(j_) = jmpFreqMat.at(i).at(j) * (1.0-occ_ave[j]);
					} else {
						freqs.at(j_) = jmpFreqMat.at(i).at(j) * (1.0-occ_eq.at(i).at(j));
					};
					freq_sum += freqs.at(j_);
					double ratio = freqs.at(j_)/freqs_prev.at(j_);
					if ( ratio < 1.0 ) {
						ratio = 1.0/ratio;
					};
					ratio_max = max(ratio_max,ratio);
				};
				dt=1.0/freq_sum;
			} while (ratio_max > 1.05);

			occ_crt[i] = 0.0;

			for (int j_ = 0 ; j_ < (adj.at(i)).size() ; j_++ ) {
				int j = adj.at(i).at(j_);
				vector<double> jmpVec_crt = {jmpVec_prev.at(0)+jmpVecMat_x.at(i).at(j),
					                           jmpVec_prev.at(1)+jmpVecMat_y.at(i).at(j),
																		 jmpVec_prev.at(2)+jmpVecMat_z.at(i).at(j)};
				if (occ_crt.count(j)) {
					occ_crt.erase(j);
				};
				nodes_crt.push_back(Node(n,j,rid,t_prev+dt,p_prev*freqs[j_]*dt,jmpVec_crt,occ_crt));
			};
			nodes_prev.pop_front();
		};

		//bool hoge = nodes_prev.empty();
		//cout << boolalpha << hoge << "\n";

		nodes_prev = nodes_crt;
		nodes_crt = {};
		
		
		// Output jmpdata.csv.modn
		jmpdata = mkjmpdata(nodes_prev);
		string n_str = boost::lexical_cast<string>(n);
		output(jmpdata, "jmpdata.csv.mod"+n_str);

		/*
		// Print out jumpdata.csv.modn
		// Combine the same jumps.
		map<string,double> jmpdata;
		for (auto itr3 = nodes_prev.begin(); itr3 != nodes_prev.end(); ++itr3) {
			int i = (*itr3).get_root_id() + 1;
			int j = (*itr3).get_site_id() + 1;
			vector<double> jmpVec = (*itr3).get_jmpVec();
			double time = (*itr3).get_time();
			double prob = (*itr3).get_prob();
			if (abs(jmpVec.at(0)) < 0.01 and abs(jmpVec.at(1)) < 0.01 and abs(jmpVec.at(2)) < 0.01) {
				;
			} else {
				stringstream ss;
				ss << fixed << i <<"," << j << ","
					 << jmpVec.at(0) << "," << jmpVec.at(1) << "," << jmpVec.at(2) << ",";
				if (jmpdata.count(ss.str())) {
					jmpdata[ss.str()] += prob/time;
				} else {
					jmpdata[ss.str()] = prob/time;
				};
			};
		};
		
		// Output files.
		string n_str = boost::lexical_cast<string>(n);
		ofstream outputfile("jmpdata.csv.mod"+n_str);
		for (auto itr4 = jmpdata.begin(); itr4 != jmpdata.end(); ++itr4) {
			string key = itr4->first;
			double val = itr4->second;
			outputfile << key << scientific << setprecision(10) << val << "\n";
		};
		outputfile.close();
		*/
	};
}


