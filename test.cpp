#include <iostream>
#include <cstdlib>
#include <string>
#include <thread>
#include "omp.h"

std::string getenv_to_string(const char *in){
	char *gotten = std::getenv(in);
	if(NULL == gotten){
		return std::string("");
	} else {
		return std::string(gotten);
	}
}

std::string getenv(const std::string& in){
	return getenv_to_string(in.c_str());
	}

int omp_thread_count(){
	int n = 0;
#pragma omp parallel reduction (+:n)
	n += 1;
	return n;
}

int main(){
	std::string envName = "OMP_NUM_THREADS";
	std::cout << envName << "    = " << getenv(envName) << std::endl;
	std::cout << "hardware_concurrency() = " << std::thread::hardware_concurrency() << std::endl; 
	std::cout << "omp_get_max_threads() = " << omp_get_max_threads() << std::endl; 
	std::cout << "omp_get_num_threads() = " << omp_get_num_threads() << std::endl;
	std::cout << "omp_thread_count() = " << omp_thread_count() << std::endl;

	//#pragma omp parallel
	#pragma omp parallel num_threads(2)
	//#pragma omp critical
	{
		std::cout << "Hello I am a thread" << omp_get_thread_num() << " of " << omp_get_num_threads() << std::endl;
		std::cout << "My C++ std::thread id is" << std::this_thread::get_id() << std::endl;
	}


	return 0;

}