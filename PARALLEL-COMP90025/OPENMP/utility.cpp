#include "utility.h"

//memory locate
int ***alloc3d(int l, int m, int n) {
    int *data = new int [l*m*n];
    int ***array = new int **[l];
    for (int i=0; i<l; i++) {
        array[i] = new int *[m];
        for (int j=0; j<m; j++) {
            array[i][j] = &(data[(i*m+j)*n]);
        }
    }
    return array;
}

//initial map, set random value
void initMap(int ***map) {
	#pragma omp parallel
	for (int i = 0; i< I; ++i){
		#pragma omp for schedule(dynamic,1)
		for (int j = 0; j < J; ++j) {
			distrubiteType(map[i][j]);
		}
	}
}


//move loop
void move(int ***map) {
	#pragma omp parallel
	for (int i = 0; i < I; ++i) {
        #pragma omp for schedule(dynamic,1)
		for (int j = 0; j < J; ++j) {

			if (map[i][j][TYPE] == EMPTY && map[i][j][MOVED] == false) {
 				map[i][j][MOVED] = true;
 			} else if (map[i][j][TYPE] != EMPTY
					   && map[i][j][MOVED] == false) {
			   seek(map, map[i][j], i, j);
			}
			
		}
	}
}

//update info after each timestamp
void updateMap(int ***map) {
    #pragma omp parallel
	for (int i = 0; i< I; ++i) {
	    #pragma omp for schedule(dynamic,1)
		for (int j = 0; j < J; ++j) {
			update(map[i][j]);
		}
	}
}

// void updateMap(int *** map) {
// #ifdef _OPENMP
// #pragma omp for schedule(dynamic,1) nowait
// #endif
// 	for (int i = 0; i< I; ++i) {
// 		for (int j = 0; j < J; ++j) {
// 			update(map[i][j]);
// 		}
// 	}
// }

// void move(int ***map) {
// #ifdef _OPENMP
// #pragma omp for schedule(dynamic,1) 
// #endif

// 	for (int ij = 0; ij< MAX; ++ij) {
// 			int x = ij/J;
// 			int y = ij%J;
// 			if (map[x][y][TYPE] == EMPTY && map[x][y][MOVED] == false) {
// 				map[x][y][MOVED] = true;
// 			} else if (map[x][y][TYPE] != EMPTY
// 					   && map[x][y][MOVED] == false) {
// 			   seek(map, map[x][y], x, y);
// 		}


// 	}
// }


//generate random number
int random(int min, int max)  {
	return min + rand() % (max - min);
}


//using to generate initial info for map randomly
void distrubiteType (int *cell) {
	int rand = random(1, 50);
//	P(rand);
    if (rand <= 4) {
		cell[TYPE]  = HUMAN;
		cell[AGE]   = random(0, HLIFETIME);
		cell[SPEED] = HSPEED;
	} else if (rand == 5) {
		cell[TYPE]  = EXPOSED;
		cell[AGE]   = random(0, ELIFETIME);
		cell[SPEED] = ESPEED;
	} else if (rand == 6) {
		cell[TYPE]  = ZOMBIE;
		cell[AGE]   = random(0, ZLIFETIME);
		cell[SPEED] = ZSPEED;
	} else {
		cell[TYPE]  = EMPTY;
		cell[AGE]   = EMPTY;
		cell[SPEED] = EMPTY;
	}
	cell[MOVED] = 0;
}

//print out result
void printResult(int ***map) {
	for (int i = 0; i< I; ++i) {
		for (int j = 0; j < J; ++j) {
			if (map[i][j][TYPE] == EMPTY)
				cout << " ";
			else 
				cout << map[i][j][TYPE];
			cout << " ";
		}
		cout << endl;
	}
}

//generate output result
void outputResult(int ***map) {
	ofstream output("output.txt");
	for (int i = 0; i< I; ++i) {
		for (int j = 0; j < J; ++j) {
			if (map[i][j][TYPE] == EMPTY)
				output << " ";
			else if  (map[i][j][TYPE] == HUMAN)
				output << "H";
			else if  (map[i][j][TYPE] == EXPOSED)
				output << "E";
			else if  (map[i][j][TYPE] == ZOMBIE)
				output << "Z";
			output << " ";
		}
		output << endl;
	}
}

//generate random coordiate
int getCoordinate(int range, int i, int size) {
	int ni = random(-(range), range);
	if (i + ni < 0)             ni = size + ni;
	else if (i + ni > size - 1) ni = i + ni - size;
	else                        ni = i + ni;
	return ni;
}

//cell seeking
void seek(int ***map,int* cell, int i,int j) {
	int ni,nj;
	int count = 0;
	if (cell[TYPE] == HUMAN) {
		while (true) {
			if (count >= 50) {
				cell[MOVED] = true;
				return;
			} else{
				ni = getCoordinate(HSPEED, i, I);
				nj = getCoordinate(HSPEED, j, J);
				if (map[ni][nj][TYPE] == EMPTY) {
					cellMove(cell, map[ni][nj]);
					return;
				}
			}
			++count;
		}
	} else if (cell[TYPE] == ZOMBIE) {
		while (true) {
			if (count >= 50) {
				cell[MOVED] = true;
				return;
			} else {
				ni = getCoordinate(ZSPEED, i, I);
				nj = getCoordinate(ZSPEED, j, J);
				if (map[ni][nj][TYPE] == EMPTY) {
					cellMove(cell, map[ni][nj]);
					return;
				} else if (map[ni][nj][TYPE]== HUMAN) {
					attack(map[ni][nj]);
					cell[MOVED] = true;
					return;
				}  
				++count;
			}
		}
	}
}

//move to new cell then reset old cell
void cellMove(int *oldC, int *newC) {
	newC[TYPE]  = oldC[TYPE];
	newC[AGE]   = oldC[AGE];
	newC[SPEED] = oldC[SPEED];
	newC[MOVED] = true;
	reset(oldC);
}

//result cell to empty
void reset(int *cell) {
	cell[TYPE]  = EMPTY;
	cell[AGE]   = EMPTY;
	cell[SPEED] = EMPTY;
}

//zombie attacks human
void attack(int *human) {
	human[TYPE]  = EXPOSED;
	human[AGE]   = 0;
	human[SPEED] = ESPEED;
	human[MOVED] = true;

}

//exposed to zombie
void exposedToZombie(int *cell) {
	cell[TYPE]  = ZOMBIE;
	cell[AGE]   = 0;
	cell[SPEED] = ZSPEED;
}

//get lifetime for different types
int getLifeTime(int *cell) {
	if (cell[TYPE] == HUMAN)
		return HLIFETIME;
	else if (cell[TYPE] == EXPOSED)
		return ELIFETIME;
	else 
		return ZLIFETIME;
}

//update cell info
void update(int *cell) {
	int isBorn;

	int lifetime = getLifeTime(cell);
	if (cell[TYPE] == EMPTY) {
		isBorn = random(0, 500);
		if (isBorn == true){
			cell[TYPE]  = HUMAN;
			cell[AGE]   = 0;
			cell[SPEED] = HSPEED;
		}
	}
	else if (cell[TYPE] == EXPOSED && ETOZ (cell[AGE]))
		exposedToZombie(cell);
	else if (DIE(cell[AGE], lifetime))
		reset(cell);
	else
		++(cell[AGE]);

	cell[MOVED] = false;
}

