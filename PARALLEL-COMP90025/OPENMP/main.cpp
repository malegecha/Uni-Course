#include "utility.h"

int main (int argc, char **argv) {

	int ***map = alloc3d(I, J, 4);
	initMap(map);
	double startTime = omp_get_wtime();
	//main loop
	for (int i = 0; i < 20; ++i){
		move(map);
		updateMap(map);
	}


	cout << "Time: " << omp_get_wtime() - startTime << "\t"
		 << omp_get_max_threads() << endl;

//	printResult(map);
//	outputResult(map);
	return 0;
}

