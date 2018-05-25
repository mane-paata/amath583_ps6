#include <iostream>
#include <cstdlib>
#include <iomanip>

#include "Vector.hpp"
#include "amath583.hpp"
#include "Timer.hpp"

int main(int argc, char* argv[]){
	if(argc !=3)
		std::cout << "Usage: " << argv[0] << "[Vector Size] [Number of times to run]" << std::endl;

	else{

		size_t dim = atoi(argv[1]);
		size_t times = atoi(argv[2]);

		double seq_time, par_time;
		
		Vector x(dim); 
		double pnorm = 0.0, snorm = 0.0;
		randomize(x);
		Timer t;

		/* sequential norm */
		t.start();
		for(int j = 0; j < times; ++j)
			snorm = twoNorm(x);
		t.stop();
		seq_time = t.elapsed()/times;

		/* parallelized norm*/
		t.start();
		for(int i = 0; i < times; ++i)
			pnorm = ompTwoNorm(x);
		t.stop();
		par_time = t.elapsed()/times;

		double speed_up = seq_time/par_time;
		double diff_norm = snorm - pnorm;
		
		/* output */
		std::cout << std::setprecision(15) << dim << "\t" << seq_time << "\t" << par_time << "\t" << speed_up << "\t" << diff_norm << std ::endl;

	}
	return 0;
}