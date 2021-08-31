#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>


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
