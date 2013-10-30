#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <cstdlib>
#ifdef _OPENMP 
#include <omp.h>
#endif
//#include <mpi.h>

using namespace std;

//constants
const int I         = 1000;
const int J         = 1500;
const int MAX       = I*J;

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
const int ELIFETIME = 10;
const int ESPEED    = 0;
const int ZLIFETIME = 200;
const int ZSPEED    = 1;


//marcos
#define P(VAR) cout << #VAR " = " << VAR << endl

#define DIE(AGE, LIFETIME) AGE == LIFETIME ? true : false    
#define ETOZ(EAGE) EAGE == ELIFETIME ? true  : false


//functions
int ***alloc3d(int l, int m, int n);
void initMap(int ***);
int random(int, int);
void distrubiteType (int *);
void printResult(int ***);
void outputResult(int ***);
int getCoordinate(int, int, int);
void reset(int *);
void seek(int ***, int *, int, int);
void move(int ***);
void cellMove(int *, int *);
void attack(int *);
void exposedToZombie(int *cell);
int getLifeTime(int *cell);
void update(int *cell);
void updateMap(int *** map);

