#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <thread>

#include "Vector.hpp"
#include "COOMatrix.hpp"
#include "amath583.hpp"

int main(int argc, char* argv[]){
	if(argc !=2){
		std::cout << "Usage: " << argv[0] << " [Size] " << std::endl;
		return -1;
	}
	size_t user_dim = atoi(argv[1]);
	size_t times = user_dim >= 100 ? 100 : 10000;
	unsigned long dim = user_dim * user_dim;
	double seq_time, par_time = 1, seq_norm, par_norm;

	Vector x(dim), y_seq(dim), y_par(dim);
	randomize(x);

	COOMatrix A(dim, dim);
	piscetize(A, user_dim, user_dim);

	/* Sequential norm */
	seq_time = seq_matvec(A, x, y_seq, times);
	seq_norm = twoNorm(y_seq);

	/* Parallel norm */
	//par_time = par_matvec(A, x, y_par, times);
	//par_norm = twoNorm(y_par);

	/* output */
	double speed_up = seq_time / par_time;
	int cores = std::thread::hardware_concurrency();
	int thread_count = 0; //omp_thread_count();
	double diff_norm = seq_norm - par_norm;
	std::cout << std::setprecision(15) << cores << "\t" << thread_count << "\t" << dim << "\t" << seq_time << "\t" << par_time << "\t" << speed_up << "\t" << diff_norm  << std ::endl;

	
	return 0;
}