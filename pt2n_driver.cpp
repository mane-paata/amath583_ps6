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
		
		Vector x(dim);
		double pnorm = 0.0, snorm = 0.0;
		randomize(x);
		Timer t1, t2;
		/////////// parallelized norm //////////
		t1.start();
		for(int i = 0; i < times; ++i)
			pnorm = ompTwoNorm(x);
		t1.stop();
		/////////// sequential norm /////////////
		t2.start();
		for(int j = 0; j < times; ++j)
			snorm = twoNorm(x);
		t2.stop();
		////////// output ////////////////////////
		std::cout << dim << "\t" << t2.elapsed()/times << "\t" << t1.elapsed()/times << "\t" << t2.elapsed()-t1.elapsed() << "\t" << snorm << "\t" << pnorm << std::setprecision(15) << std ::endl;

	}
	return 0;
}