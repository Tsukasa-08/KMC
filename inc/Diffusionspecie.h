#include <vector>
#include <cmath>


//拡散種クラスを定義
class Diffusionspecie {

private:
	int diffusion_id;
	std::vector<double> jump_total;

public:

	Diffusionspecie() : diffusion_id(-1), jump_total(3,0.0) {
	};

	//セッタ
	void set_diffusion_id(int d_id) { diffusion_id = d_id; }

	void set_jump_total(std::vector<double> jt) { jump_total = jt; }

	//ゲッタ
	int get_diffusion_id() { return diffusion_id; }

	std::vector<double> get_jump_total() { return jump_total; }

	//拡散係数を算出する関数(引数にtotal_time)
	std::vector<double> get_D(double t) {
		std::vector<double> D_3d(3,0.0);
		for (int i = 0; i != D_3d.size(); i++) { 
			double mean_square_displacement = pow(get_jump_total()[i],2.0);
			D_3d[i] = mean_square_displacement/(2*t);
		}
		return D_3d;
	}

};
