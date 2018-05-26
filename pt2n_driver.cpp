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

    double seq_time = 0.0, par_time = 0.0;
    
    Vector x(dim); 
    double pnorm = 0.0, snorm = 0.0;
    randomize(x);
    Timer t;

    /* sequential norm */
    for(int j = 0; j < times; ++j)
		{
		  t.start();
      snorm = twoNorm(x);
      t.stop();
      seq_time += t.elapsed();
		}
		seq_time = seq_time / times;

    /* parallelized norm*/
    for(int j = 0; j < times; ++j)
		{
		  t.start();
      pnorm = ompTwoNorm(x);
      t.stop();
      par_time += t.elapsed();
		}
		par_time = par_time / times;

    double speed_up = seq_time/par_time;
    double diff_norm = std::abs(snorm - pnorm);
    
    /* output */
    std::cout << std::setprecision(15) << dim << "\t" << seq_time << "\t" << par_time << "\t" << speed_up << "\t" << diff_norm << std ::endl;

  }
  return 0;
}