#include <iostream>

#include "Vector.hpp"
#include "COOMatrix.hpp"
#include "amath583.hpp"

int main(int argc, char* argv[]){
	if(argc !=2){
		std::cout << "Usage: " << argv[0] << " [Size] " << std::endl;
		return -1;
	}
	size_t user_dim = atoi(argv[1]);
	driver_helper(user_dim, true);
	return 0;
}