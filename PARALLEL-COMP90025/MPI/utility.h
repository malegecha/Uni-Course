#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <assert.h>
#ifdef _OPENMP 
#include <omp.h>
#endif
#include <mpi.h>
#include <stdio.h>
#include "Image.h"
using namespace std;

//constants
const int LEFT      = 0;
const int RIGHT     = 1;
const int UP        = 2;
const int DOWN      = 3;

const int MAPI      = 1000;
const int MAPJ      = 1500;
const int TIME      = 200;
const int BIRTHRATE = 365*80;

const int TYPE      = 0;
const int AGE       = 1;
const int SPEED     = 2;
const int MOVED     = 3;

const int EMPTY     = -1;
const int HUMAN     = 1;
const int EXPOSED   = 2;
const int ZOMBIE    = 3;

const int HLIFETIME = 80*365;  
const int HSPEED    = 1;   
const int ELIFETIME = 50;
const int ESPEED    = 0;
const int ZLIFETIME = 200;
const int ZSPEED    = 1;


//marcos
#define P(VAR) cout << #VAR " = " << VAR << endl
#define DIE(AGE, LIFETIME) AGE == LIFETIME ? true : false    
#define ETOZ(EAGE) EAGE == ELIFETIME ? true  : false

//functions
int ***alloc3d(int, int, int);
int ****alloc4d(int, int, int, int);
void initMap(int ***, int, int);
int random(int, int);
void distrubiteType (int *);
int *getCoordinate(int, int, int);
void reset(int *);
void seek(int ***, int *, int, int);
void move(int ***, int, int);
void cellMove(int *, int *);
void attack(int *);
void exposedToZombie(int *);
void update(int *);
void updateMap(int ***, int, int);
void freeMap4d(int ****, int, int);
void freeMap3d(int ***, int, int);
void decisionForBoundary(int ***, int ***, int, int, int, int);
void calPup(int ****, int, int);
int getRank();
int getL(int);
int getR(int);
int getT(int);
int getB(int);
int getIntArg(int argc, char *argv[], const char *arg, int defaultValue);
int getSize();
int getLifeTime(int *);
void writeTGA(int ***, int ***, int ***, int ***, int , int, int);
