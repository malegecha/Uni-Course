#include "utility.h"
//memory locate
int ***alloc3d(int l, int m, int n) {
	int *data = new int[l * m * n];
	int ***array = new int **[l];
	for (int i = 0; i < l; i++) {
		array[i] = new int *[m];
		for (int j = 0; j < m; j++) {
			array[i][j] = &(data[(i * m + j) * n]);
		}
	}
	return array;
}

int ****alloc4d(int D1, int D2,int D3,int D4) {
	int**** v;
	v = new int***[D1];
	for (int i = 0; i< D1; ++i)
		v[i] = alloc3d(D2, D3, D4);
	return v;

}
void initMap(int ***map, int I, int J) {
	for (int i = 1; i < I + 1; ++i) {
		for (int j = 1; j < J + 1; ++j) {
			distrubiteType(map[i][j]);
		}
	}
 }

//move loop
void move(int ***map, int I, int J) {
	for (int i = 1; i < I + 1; ++i) {
		for (int j = 1; j < J + 1; ++j) {
			if (map[i][j][TYPE] == EMPTY && map[i][j][MOVED] == false) {
				map[i][j][MOVED] = true;
			} else if (map[i][j][TYPE] != EMPTY && map[i][j][MOVED] == false) {
				seek(map, map[i][j], i, j);
			}
		}
	}
}

//update info after each timestamp
void updateMap(int ***map, int I, int J) {
	for (int i = 1; i < I + 1; ++i) {
		for (int j = 1; j < J + 1; ++j) {
			update(map[i][j]);
		}
	}
}

//generate random number
int random(int min, int max) {
	return min + rand() % (max - min);
}

//using to generate initial info for map randomly
void distrubiteType(int *cell) {
	int rand = random(1, 1000);
	if (rand <= 34) {
		cell[TYPE] = HUMAN;
		cell[AGE] = random(0, HLIFETIME);
		cell[SPEED] = HSPEED;
	}  else if (rand == 36 || rand == 35) {
		cell[TYPE] = ZOMBIE;
		cell[AGE] = random(0, ZLIFETIME);
		cell[SPEED] = ZSPEED;
	} else {
		cell[TYPE] = EMPTY;
		cell[AGE] = EMPTY;
		cell[SPEED] = EMPTY;
	}
	cell[MOVED] = 0;
}


int *getCoordinate(int range, int i, int j) {
	int r = random(1, 10);
	int *xy = new (int);
	if (r == 1) {
		xy[0] = i + range;
		xy[1] = j;
	} else if (r == 2) {
		xy[0] = i - range;
		xy[1] = j;
	} else if (r == 3) {
		xy[0] = i ;
		xy[1] = j + range;
	} else if (r == 4) {
		xy[0] = i;
		xy[1] = j - range;
	} else {
		xy[0] = i;
		xy[1] = j;
	}
	return xy;
}
//cell seeking
void seek(int ***map, int* cell, int i, int j) {
	int *xy;
	int count = 0;
	if (cell[TYPE] == HUMAN) {
		while (true) {
			if (count >= 4) {
				cell[MOVED] = true;
				return;
			} else {
				xy = getCoordinate(HSPEED, i, j);
				if (map[xy[0]][xy[1]][TYPE] == EMPTY) {
					cellMove(cell, map[xy[0]][xy[1]]);
					return;
				} else if (map[xy[0]][xy[1]][TYPE] == HUMAN) {
//					cellMove(cell, map[xy[0]][xy[1]]);
					xy = getCoordinate(HSPEED, i, j);
					if (map[xy[0]][xy[1]][TYPE] == EMPTY) {
						int isBorn = random(0, BIRTHRATE);
						if (isBorn == true) {
							map[xy[0]][xy[1]][TYPE]	 = HUMAN;
							map[xy[0]][xy[1]][AGE] = 0;
							map[xy[0]][xy[1]][SPEED] = HSPEED;
							map[xy[0]][xy[1]][MOVED] = true;
						}
						cell[MOVED] = true;
					}
					return;
				}
			}
			++count;
		}
	} else if (cell[TYPE] == ZOMBIE) {
		while (true) {
			if (count >= 4) {
				cell[MOVED] = true;
				delete xy;
				return;
			} else {
				xy = getCoordinate(ZSPEED, i, j);
				if (map[xy[0]][xy[1]][TYPE] == EMPTY) {
					cellMove(cell, map[xy[0]][xy[1]]);
					return;
				} else if (map[xy[0]][xy[1]][TYPE] == HUMAN) {
					attack(map[xy[0]][xy[1]]);
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
	newC[TYPE] = oldC[TYPE];
	newC[AGE] = oldC[AGE];
	newC[SPEED] = oldC[SPEED];
	newC[MOVED] = true;
	reset(oldC);
}

//result cell to empty
void reset(int *cell) {
	cell[TYPE] = EMPTY;
	cell[AGE] = EMPTY;
	cell[SPEED] = EMPTY;
}

//zombie attacks human
void attack(int *human) {
	human[TYPE] = EXPOSED;
	human[AGE] = 0;
	human[SPEED] = ESPEED;
}

//exposed to zombie
void exposedToZombie(int *cell) {
	cell[TYPE] = ZOMBIE;
	cell[AGE] = 0;
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
	int lifetime = getLifeTime(cell);

	if (cell[TYPE] == EXPOSED && ETOZ(cell[AGE]))
		exposedToZombie(cell);
	else if (DIE(cell[AGE], lifetime))
		reset(cell);
	else
		++(cell[AGE]);
	cell[MOVED] = false;
}

void freeMap4d(int ****map, int I, int J) {
	int size = getSize();
	for (int i = 0; i < size; ++i)
		freeMap3d(map[i], I, J);
}
void freeMap3d(int ***map, int I, int J) {
	for (int i = 0; i < I + 2; ++i) {
		for (int j = 0; j < J + 2; ++j)
			delete[] map[i][j];

		delete[] map[i];
	}
	delete[] map;
}


void calPup(int ****map, int I, int J) {
	int num_human, num_exposed, num_zombie;
	num_human = num_exposed = num_zombie = 0;
	ofstream out;
	out.open("stat.csv", ios::app);
	int size = getSize();
	for (int k = 0; k < size; ++k){
		for (int i = 1; i < I + 1; ++i) {
			for (int j = 1; j < J + 1; ++j) {
				if (map[k][i][j][TYPE] == HUMAN){
					num_human++;
				} else if (map[k][i][j][TYPE] == EXPOSED) {
					num_exposed++;
				} else if (map[k][i][j][TYPE] == ZOMBIE) {
					num_zombie++;
				} 
			}
		}
	}
	out << num_exposed<< ",\t" <<num_human << ",\t" << num_zombie  <<"\n";
	out.close();
}


//boundary strategy
void decisionForBoundary(int ***buffer, int ***map, int type, int k, int I, int J) {
	if (type == I) {

		for (int i = 0; i < I ; i++) {
			if(buffer[i][0][MOVED] == true && buffer[i][0][0]!=EMPTY && map[i+1][k][0] == EMPTY) {
				map[i+1][k][0] = buffer[i][0][0];
				map[i+1][k][1] = buffer[i][0][1];
				map[i+1][k][2] = buffer[i][0][2];
				map[i+1][k][3] = buffer[i][0][3];
			} else if (buffer[i][0][0] == map[i+1][k][0]){
				map[i+1][k][0] = buffer[i][0][0];
				map[i+1][k][1] = buffer[i][0][1];
				map[i+1][k][2] = buffer[i][0][2];
				map[i+1][k][3] = buffer[i][0][3];
			}
		}
	}else {
		for (int i = 0; i < J ; i++) {
			if(buffer[0][i][MOVED] == true && buffer[0][i][0]!=EMPTY && map[k][i+1][0] ==EMPTY) {
				map[k][i+1][0] = buffer[0][i][0];
				map[k][i+1][1] = buffer[0][i][1];
				map[k][i+1][2] = buffer[0][i][2];
				map[k][i+1][3] = buffer[0][i][3];
			} else if (buffer[0][i][0] == map[k][i+1][0]){
				map[k][i+1][0] = buffer[0][i][0];
				map[k][i+1][1] = buffer[0][i][1];
				map[k][i+1][2] = buffer[0][i][2];
				map[k][i+1][3] = buffer[0][i][3];
			}
		}
	}
}

// get input arguments
int getIntArg(int argc, char *argv[], const char *arg, int defaultValue) {
    for (int i = 1; i < argc - 1; i ++) {
        if (strcmp(argv[i], arg) == 0) {
            return(atoi(argv[i+1]));
        }
    }
    return(defaultValue);
}

//get MPI size
int getSize() {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return(size);
}

//get MPI rank
int getRank() {
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	return world_rank;
}

//get left neighbour
int getL(int j) {
	int rank = getRank();
	if (rank  % j == 0)
		return rank + j -1;
	return rank - 1;
}

//get right neighbour
int getR(int j) {
	int rank = getRank();
	if ((rank + 1) % j == 0)
		return rank - j +1;
	return rank  +1;
}

//get top neighbour
int getT(int j) {
	int rank = getRank();
	if ((rank - j) < 0)
		return rank - j + getSize();
	return rank - j;
}

//get bottom neighbour
int getB(int j) {
	int rank = getRank();
	if((rank + j) > getSize() - 1)
		return rank + j - getSize();
	return rank + j;
}

void writeTGA(int ***map0, int ***map1, int ***map2, int ***map3, int fileNumber, int I, int J) {

	TGAImage *tgaimage = new TGAImage(I*2, J*2);
	Colour human;
	human.r = 0;
	human.g = 255;
	human.b = 0;
	human.a = 255;
	
	Colour empty;
	empty.r = 0;
	empty.g = 0;
	empty.b = 0;
	empty.a = 255;

	Colour exposed;
	exposed.r = 0;
	exposed.g = 0;
	exposed.b = 255;
	exposed.a = 255;

	Colour zombie;
	zombie.r = 255;
	zombie.g = 0;
	zombie.b = 0;
	zombie.a = 255;

	
	char filename[30];
	int num_human, num_exposed, num_zombie;
	num_human = num_exposed = num_zombie = 0;
	ofstream out;
	out.open("stat.csv", ios::app);
	for (int i = 1; i < I + 1; ++i) {
		for (int j = 1; j < J + 1; ++j) {
			if (map0[i][j][TYPE] == HUMAN){
				tgaimage->setPixel(human, i-1, j-1);
				num_human++;
			} else if (map0[i][j][TYPE] == EXPOSED) {
				tgaimage->setPixel(exposed, i-1, j-1);
				num_exposed++;
			} else if (map0[i][j][TYPE] == ZOMBIE) {
				tgaimage->setPixel(zombie, i-1, j-1);
				num_zombie++;
			}  else if (map0[i][j][TYPE] == EMPTY) {
				tgaimage->setPixel(empty, i-1, j-1);
			}
		}
		for (int j = 1; j < J + 1; ++j) {
			if (map1[i][j][TYPE] == HUMAN) {
				tgaimage->setPixel(human, i-1, J+j-1);
				num_human++;
			} else if (map1[i][j][TYPE] == EXPOSED) {
				tgaimage->setPixel(exposed, i-1, J+j-1);
				num_exposed++;
			} else if (map1[i][j][TYPE] == ZOMBIE){
				tgaimage->setPixel(zombie, i-1, J+j-1);
				num_zombie++;
			} else if (map1[i][j][TYPE] == EMPTY) {
				tgaimage->setPixel(empty, i-1, J+j-1);
			}
		}
	}
	for (int i = 1; i < I+1; ++i) {
		for (int j = 1; j < J + 1; ++j) {
			if (map2[i][j][TYPE] == HUMAN){
				tgaimage->setPixel(human, I+i-1, j-1);
				num_human++;
			} else if (map2[i][j][TYPE] == EXPOSED) {
				tgaimage->setPixel(exposed, I+i-1, j-1);
				num_exposed++;
			} else if (map2[i][j][TYPE] == ZOMBIE) {
				tgaimage->setPixel(zombie, I+i-1, j-1);
				num_zombie++;
			} else if (map2[i][j][TYPE] == EMPTY) {
				tgaimage->setPixel(empty, I+i-1, j-1);
			}
		}
		for (int j = 1; j < J + 1; ++j) {
			if (map3[i][j][TYPE] == HUMAN) {
				tgaimage->setPixel(human, I+i-1, J+j-1);
				num_human++;
			} else if (map3[i][j][TYPE] == EXPOSED){
				tgaimage->setPixel(exposed, I+i-1, J+j-1);
				num_exposed++;
			} else if (map3[i][j][TYPE] == ZOMBIE) {
				tgaimage->setPixel(zombie, I+i-1, J+j-1);
				num_zombie++;
			} else if (map3[i][j][TYPE] == EMPTY) {
				tgaimage->setPixel(empty, I+i-1, J+j-1);
			}
		}
	}
	out << num_exposed<< ",\t" <<num_human << ",\t" << num_zombie  <<"\n";
	
	sprintf(filename, "%d%s", fileNumber, ".tga");
	tgaimage->WriteImage(filename);
	out.close();
}
