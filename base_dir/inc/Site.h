#include <vector>
#include "Jump.h"

//サイトクラスを定義
class Site {

private:
	int site_id;
	int site_atom;
	int diffusion_id;
	std::vector<double> site_frac_coords;
	std::vector<Jump> jumps_from_here;
	std::vector<Jump> jumps_to_here;
	std::set<int> blocking_mate_list;

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

	void set_blocking_mate_list(std::set<int> st) { blocking_mate_list = st ; } ;

	//ゲッタ
	int get_site_id() { return site_id; } ;

	int get_site_atom() { return site_atom; } ;

	int get_diffusion_id() { return diffusion_id; } ;

	std::vector<double> get_site_frac_coords() { return site_frac_coords; } ;

	std::vector<Jump> get_jumps_from_here() { return jumps_from_here; } ;

	std::vector<Jump> get_jumps_to_here() { return jumps_to_here; } ;

	std::set<int> get_blocking_mate_list() { return blocking_mate_list; } ;
};
