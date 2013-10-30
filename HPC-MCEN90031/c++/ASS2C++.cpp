/////////////////////////////////////////////////////////////////////////////
//
//  MCEN90031 Applied High Performance Computing Assignment 2
//
//  Problem:		    The Heat Equation
//
//  Compilation:		g++ ASS2C++.cpp -o ASS2C++
//
//  Execution:			./ASS2C++ Box.grid
//
/////////////////////////////////////////////////////////////////////////////


#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <cstring>
#include <iomanip>
using namespace std;

// Class definitions
class SparseMatrix {
public:
	SparseMatrix(int M, int nnz) {
		nnz_		= 0;
		M_			= M;
		allocSize_	= nnz;
		val_		= new double	[allocSize_];
		col_		= new int		[allocSize_];
		row_		= new int		[M_+1];
		
		memset(val_, 0.0, allocSize_*sizeof(double));
		memset(col_, 0,   allocSize_*sizeof(int));
		memset(row_, 0,	  (M_+1)    *sizeof(int));
	}
	SparseMatrix() {
		val_		= NULL;
		col_		= NULL;
		row_		= NULL;
		nnz_		= 0;
		M_			= 0;
		allocSize_	= 100;
	}
	~SparseMatrix() {
		if(val_)	delete	[] val_;
		if(col_)	delete	[] col_;
		if(row_)	delete	[] row_;
	}
	void initialize(int M, int nnz) {
		nnz_		= 0;
		M_			= M;
		allocSize_	= nnz;
		val_		= new double	[allocSize_];
		col_		= new int		[allocSize_];
		row_		= new int		[M_+1];
		
		memset(val_, 0.0, allocSize_*sizeof(double));
		memset(col_, 0,   allocSize_*sizeof(int));
		memset(row_, 0,	  (M_+1)    *sizeof(int));
		return;
	}
	inline double& operator()(int m, int n) {
		int		k			= row_[m];
		bool	foundCol	= false;
		int		insertPos	= k;
		
		while(k<row_[m+1] && !foundCol) {
			if(col_[k]==n)
				foundCol=true;
			if(col_[k]< n)
				insertPos++;
			k++;
		}
		if(foundCol) {
			return val_[k-1];
		} else {
			if(nnz_+1>=allocSize_)
				this->reallocate();
			for(int i=nnz_; i>insertPos; i--) {
				col_[i] = col_[i-1];
				val_[i] = val_[i-1];
			}
			
			nnz_++;
			for(int i=m+1; i<M_; i++)
				row_[i]++;
			row_[M_]		= nnz_;
			col_[insertPos]	= n;
			val_[insertPos]	= 0.0;
			
			return val_[insertPos];
		}
	}
	inline double& operator()(int k) {
			return val_[k];
	}
	void operator = (const SparseMatrix& A) {
			nnz_		= A.nnz_;
			M_			= A.M_;
			allocSize_	= A.allocSize_;
			val_		= new double	[allocSize_];
			col_		= new int		[allocSize_];
			row_		= new int		[M_+1];
		
			memcpy(val_, A.val_, allocSize_*sizeof(double));
			memcpy(col_, A.col_, allocSize_*sizeof(int));
			memcpy(row_, A.row_, (M_+1)    *sizeof(int));
	}
	inline void	multiply(double* u, double* v, int M=0) {
			if(M==0)
				M=M_;
			for(int m=0; m<M; m++) {
				u[m] = 0.0;
				for(int k=row_[m]; k<row_[m+1]; k++) {
					if(col_[k]<=M) 
						u[m] += val_[k]*v[col_[k]];
				}
			}
			return;
		}
	void print() {
			for(int k=0; k<nnz_; k++)
				cout << val_[k] << endl;
			return;
	}
	void write() {
		fstream file;
		file.open("A.data", ios::out);
		for(int k=0; k<nnz_; k++)
			file << val_[k] << endl;
		for(int k=0; k<nnz_; k++)
			file << col_[k]+1 <<endl;
		for(int m=0; m<M_; m++)
			file << row_[m]+1 << endl;
		file.close();
		return;
	}
	int	getNNZ() {
		return nnz_;
	}
protected:
	void reallocate() {		
		cout << "\nReallocating!!!\n";
		// Double the memory allocation size
		allocSize_ += allocSize_;
		
		// Create some temporary arrays of the new size
		double* tempVal = new double	[allocSize_];
		int*	tempCol = new int		[allocSize_];
		
		// Copy all of the data from the old to the new arrays
		memcpy(tempVal, val_, nnz_*sizeof(double));
		memcpy(tempCol, col_, nnz_*sizeof(int));
		
		// Delete the memory used by the old arrays
		delete [] val_;
		delete [] col_;
		
		// Assign the addresses of the new arrays
		val_	= tempVal;
		col_	= tempCol;
		
			return;
		}
private:
	double*	val_;	// [nnz] Stores the nonzero elements of the matrix
	int*	col_;	// [nnz] Stores the column indices of the elements in each row
	int*	row_;	// [M+1] Stores the locations in val that start a row
	int		nnz_;
	int		M_;
	int		allocSize_;
};

class	Boundary {
public:
	Boundary() {
	}
	string  name_;
    string  type_;
	int		N_;
    int*    indices_;
	double	value_;
};

// Function declarations
void	getGrid(char* filename);
void	assembleSystem();
void	solveWithConjugateGradient(SparseMatrix& A, double* phi, double* b, int M);
void writeData(int l);	
// Global variables
double			t_min		=    0;
double			t_max		=  100;
double			Delta_t		=    1;
const double	rho         = 8954;			
const double	C			=  380;			
const double	k			=  386;			
const double	h			=  100;			
const double	T_air		=  300;			
const double	Q_cpu		= 40000;		
const int		N_t			= (t_max-t_min)/Delta_t+1;
int				N_p			= 0;
int				N_f			= 0;
int				N_e			= 0;
int				N_b			= 0;
double**		Points		= NULL;
int**			Faces		= NULL;
int**			Elements	= NULL;
Boundary*		Boundaries	= NULL;
double*			Omega		= NULL;
double*			Gamma		= NULL;
double*			T			= NULL;
SparseMatrix	M;
SparseMatrix	K;
double*			s			= NULL;
SparseMatrix	A;
double*			b			= NULL;

int	main(int argc, char** argv) {
	getGrid(argv[1]);
	// Allocate arrays
	T			= new double[N_p];
	s			= new double[N_p];
	b			= new double[N_p];
	
	// Set initial condition
	for(int m=0; m<N_p; m++)
		T[m]	= T_air;
	M.initialize(N_p, 20*N_p);
	K.initialize(N_p, 20*N_p);
	memset(s, 0.0, N_p*sizeof(double));
	assembleSystem();
	A		= M;
	for(int k=0; k<A.getNNZ(); k++)
		A(k) -= Delta_t*(K(k));
			
	fstream file;
	file.open("ass2.data", ios::out);
	
	// Main time marching loop
	for(int l=0; l<N_t-1; l++) {
		cout << "l = " << l;
		// Assemble b
		M.multiply(b, T);
		for(int m=0; m<N_p; m++)
			b[m]	+= Delta_t*s[m];

		solveWithConjugateGradient(A, T, b, N_p);
		writeData(l);
	}
		
	file.close();
	// Delete arrays
	for(int p=0; p<N_p; p++)
		delete [] Points[p];
	for(int f=0; f<N_f; f++)
		delete [] Faces[f];
	for(int e=0; e<N_e; e++)
		delete [] Elements[e];
	delete [] Points;
	delete [] Faces;
	delete [] Elements;
	delete [] Boundaries;
	delete [] Gamma;
	delete [] Omega;
	delete [] T;
	delete [] s;
	delete [] b;
	
	return 0;
}

void getGrid(char* filename) {
	fstream		file;
	char		tempName[64];
	cout << "Reading input file...\n";
	
	file.open(filename);
	file >> tempName >> N_p;
	file >> tempName >> N_f;
	file >> tempName >> N_e;
	file >> tempName >> N_b;
	
	Points			= new double*	[N_p];
	Faces			= new int*		[N_f];
	Elements		= new int*		[N_e];
	Boundaries		= new Boundary  [N_b];
	
	file >> tempName;

	for(int p=0; p<N_p; p++) {
		Points[p] = new double [3];
		file >> Points[p][0] >> Points[p][1] >> Points[p][2];
	}
	
	file >> tempName;
	for(int f=0; f<N_f; f++) {
		Faces[f] = new int [3];
		file >> Faces[f][0] >> Faces[f][1] >> Faces[f][2];
	}
	file >> tempName;
	for(int e=0; e<N_e; e++) {
		Elements[e] = new int [4];
		file >> Elements[e][0] >> Elements[e][1]
			 >> Elements[e][2] >> Elements[e][3];
	}
	file >> tempName;

	for(int b=0; b<N_b; b++) {
		file >> Boundaries[b].name_ >> Boundaries[b].type_ >> Boundaries[b].N_;
        Boundaries[b].indices_  = new int [Boundaries[b].N_];
        for(int n=0; n<Boundaries[b].N_; n++)
            file >> Boundaries[b].indices_[n];
        file >> Boundaries[b].value_;
	}
	file.close();
	return;
}

void assembleSystem() {	
	double	x[4];
	double	y[4];
	double  z[4];
	double	G[3][4];
	double	Gp[3]		= {0.0, 0.0, 0.0};
	double	Gq[3]		= {0.0, 0.0, 0.0};
	double	M_e[4][4]	= {{2.0, 1.0, 1.0, 1.0},
						   {1.0, 2.0, 1.0, 1.0},
						   {1.0, 1.0, 2.0, 1.0},
						   {1.0, 1.0, 1.0, 2.0}};

	double	K_f[3][3]	= {{2.0, 1.0, 1.0},
						   {1.0, 2.0, 1.0},
						   {1.0, 1.0, 2.0}};
	double	s_e[3]		= {1.0, 1.0, 1.0};
	int		Nodes[4]	= {0, 0, 0, 0};
	int		m;
	int		n;
	
	cout << "Assembling system..." << flush;
	
	Gamma       = new double	[N_f];
	Omega		= new double	[N_e];

    // Calculate face lengths
   	cout << "Calculate face lengths... " << endl;
	// Calculate face lengths
    for(int f=0; f<N_f; f++) {
		for(int p=0; p<3; p++) {
			x[p]	= Points[Faces[f][p]][0];
			y[p]	= Points[Faces[f][p]][1];
			z[p]	= Points[Faces[f][p]][2];
		}
		Gamma[f] = (sqrt( pow(((y[1]-y[0])*(z[2]-z[0]) - (z[1]-z[0])*(y[2]-y[0])),2) 
						  +pow(((z[1]-z[0])*(x[2]-x[0]) - (x[1]-x[0])*(z[2]-z[0])),2)
						  +pow(((x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0])),2)))/2;
	}

	cout << "Calculate element areas... " << endl;
    // Calculate element areas
    for(int e=0; e<N_e; e++) {
		for(int p=0; p<4; p++) {
			x[p]	= Points[Elements[e][p]][0];
			y[p]	= Points[Elements[e][p]][1];
			z[p]	= Points[Elements[e][p]][2];
		}
		//	Omega[e]    = ( (x[0]*y[1]-x[1]*y[0]) + (x[2]*y[0])-(x[0]*y[2]) + (x[1]*y[2]-x[2]*y[1]) )/2;
		Omega[e] = (abs(   x[0]*y[1]*z[2] - x[0]*y[2]*z[1] - x[1]*y[0]*z[2] 
						   + x[1]*y[2]*z[0] + x[2]*y[0]*z[1] - x[2]*y[1]*z[0] 
						   - x[0]*y[1]*z[3] + x[0]*y[3]*z[1] + x[1]*y[0]*z[3] 
						   - x[1]*y[3]*z[0] - x[3]*y[0]*z[1] + x[3]*y[1]*z[0] 
						   + x[0]*y[2]*z[3] - x[0]*y[3]*z[2] - x[2]*y[0]*z[3] 
						   + x[2]*y[3]*z[0] + x[3]*y[0]*z[2] - x[3]*y[2]*z[0] 
						   - x[1]*y[2]*z[3] + x[1]*y[3]*z[2] + x[2]*y[1]*z[3] 
						   - x[2]*y[3]*z[1] - x[3]*y[1]*z[2] + x[3]*y[2]*z[1] )) /6;
	}
	
	cout << "Assemble M, K, and s... " << endl;
	// Assemble M, K, and s
	for(int e=0; e<N_e; e++) {
		for(int p=0; p<4; p++)
			Nodes[p] = Elements[e][p];
		for(int p=0; p<4; p++) {
			x[p]	= Points[Nodes[p]][0];
			y[p]	= Points[Nodes[p]][1];
			z[p]	= Points[Nodes[p]][2];
		}
		
		G[0][0]	= (y[3]-y[1])*(z[2]-z[1]) - (y[2]-y[1])*(z[3]-z[1]);
		G[0][1]	= (y[2]-y[0])*(z[3]-z[2]) - (y[2]-y[3])*(z[0]-z[2]);
		G[0][2]	= (y[1]-y[3])*(z[0]-z[3]) - (y[0]-y[3])*(z[1]-z[3]);
		G[0][3] = (y[0]-y[2])*(z[1]-z[0]) - (y[0]-y[1])*(z[2]-z[0]);
		G[1][0]	= (x[2]-x[1])*(z[3]-z[1]) - (x[3]-x[1])*(z[2]-z[1]);
		G[1][1]	= (x[3]-x[2])*(z[2]-z[0]) - (x[0]-x[2])*(z[2]-z[3]);
		G[1][2]	= (x[0]-x[3])*(z[1]-z[3]) - (x[1]-x[3])*(z[0]-z[3]);
		G[1][3] = (x[1]-x[0])*(z[0]-z[2]) - (x[2]-x[0])*(z[0]-z[1]);
		G[2][0]	= (x[3]-x[1])*(y[2]-y[1]) - (x[2]-x[1])*(y[3]-y[1]);
		G[2][1]	= (x[2]-x[0])*(y[3]-y[2]) - (x[2]-x[3])*(y[0]-y[2]);
		G[2][2]	= (x[1]-x[3])*(y[0]-y[3]) - (x[0]-x[3])*(y[1]-y[3]);
		G[2][3] = (x[0]-x[2])*(y[1]-y[0]) - (x[0]-x[1])*(y[2]-y[0]);

		// Outer loop over each node
		for(int p=0; p<4; p++) {
			m		= Nodes[p];
			Gp[0]	= G[0][p];
			Gp[1]	= G[1][p];
			Gp[2]   = G[2][p];
			// Inner loop over each node
			for(int q=0; q<4; q++) {
				n			= Nodes[q];
				Gq[0]		= G[0][q];
				Gq[1]		= G[1][q];
				Gq[2]       = G[2][q];
				M(m,n)		+=	(M_e[p][q]*rho*C*Omega[e]/20);
				
				K(m,n)		-=	(k*(Gp[0]*Gq[0]+Gp[1]*Gq[1]+Gp[2]*Gq[2])/(36*Omega[e]));
			}
		}
	}
	cout << "Apply boundary conditions... " << flush;
	// Apply boundary conditions
    for(int b=0; b<N_b; b++) {
		if (Boundaries[b].type_=="neumann") {
			for(int f=0; f<Boundaries[b].N_; f++) {
				for(int p=0; p<3; p++) {
					Nodes[p]	= Faces[Boundaries[b].indices_[f]][p];
					m			= Nodes[p];
					s[m]		+=s_e[p]*Q_cpu*Gamma[Boundaries[b].indices_[f]]/3;
				}
			}
		} else if (Boundaries[b].type_=="robin") {
			for(int f=0; f<Boundaries[b].N_; f++) {
				for(int p=0; p<3; p++) {
					Nodes[p]    = Faces[Boundaries[b].indices_[f]][p];
					m			= Nodes[p];
					for(int q=0; q<3; q++) {
						Nodes[q]= Faces[Boundaries[b].indices_[f]][q];
						n		= Nodes[q];
						K(m,n)	-= K_f[p][q]*h*Gamma[Boundaries[b].indices_[f]]/12;
					}
					s[m]		+= s_e[p]*h*T_air*Gamma[Boundaries[b].indices_[f]]/3;
				}
			}
		}
	}
	cout << "Done.\n" << flush;
	return;
}

void	solveWithConjugateGradient(SparseMatrix& A, double* phi, double* b, int M) {
	double*	r0				= new double [M];
	double*	rn				= new double [M];
	double*	d				= new double [M];
	double*	Ad				= new double [M];
	double*	Aphi			= new double [M];
	double	alpha			= 0.0;
	double	beta			= 0.0;
	double	r_norm			= 0.0;
	double	tolerance		= 1e-8;
	double	maxIterations	= 1e+3;
	double	r0Tr0			= 0.0;
	double	rnTrn			= 0.0;
	double	dTAd			= 0.0;
	int		k				= 0;
	int		m				= 0;
	int		n				= 0;
	int		N				= M;
	
	memset(r0, 0.0, M*sizeof(double));
	memset(rn, 0.0, M*sizeof(double));
	memset(d,  0.0, M*sizeof(double));
	memset(Ad, 0.0, M*sizeof(double));
	
	// Compute the initial residual and it's norm
	A.multiply(Aphi, phi, M);
	for(m=0; m<M; m++)
	{
		r0[m]	= b[m] - Aphi[m];
		d[m]	= r0[m];
		r0Tr0  += r0[m]*r0[m];
	}
	r_norm = sqrt(r0Tr0);
	
	// Main iterative loop
	while(r_norm>tolerance && k<maxIterations)
	{
		dTAd	= 0.0;
		A.multiply(Ad, d, M);
		for(m=0; m<M; m++)
		{
			dTAd   += d[m]*Ad[m];
		}
		alpha  	= r0Tr0/dTAd;
		for(m=0; m<M; m++)
		{
			phi[m] += alpha*d[m];
		}
		for(m=0; m<M; m++)
		{
			rn[m]	= r0[m] - alpha*Ad[m];
		}
		rnTrn	= 0.0;
		for(m=0; m<M; m++)
		{
			rnTrn  += rn[m]*rn[m];
		}
		beta  	= rnTrn/r0Tr0;
		for(m=0; m<M; m++)
		{
			d[m] = rn[m] + beta*d[m];
		}
		for(m=0; m<M; m++)
		{
			r0[m] = rn[m];
		}
		r0Tr0	= rnTrn;
		r_norm	= sqrt(rnTrn);
		k++;   
	}
	
	cout << ", k = " << k << ", r_norm = " << r_norm << endl;
	
	delete [] r0;
	delete [] rn;
	delete [] d;
	delete [] Ad;
	delete [] Aphi;
	
	return;
}

	
void writeData(int l) {	
	fstream         file;
    char            fileName[64];
    
    sprintf(fileName, "heat_%04d.vtk", l);
	
    file.open(fileName, ios::out);
	
    file << "# vtk DataFile Version 2.0"<< endl;
    file << "untitled, Created by me"	<< endl;
    file << "ASCII"						<< endl;
    file << "DATASET UNSTRUCTURED_GRID"	<< endl;
    
    file << "POINTS " << N_p << " double" << endl;
    for(int p=0; p<N_p; p++)
	file << setw(6) << setprecision(5) << fixed << Points[p][0] << "\t" << Points[p][1] << "\t" << Points[p][2] << endl;
    
    file << "CELLS " << N_e << " " << 5*N_e << endl;
    for(int e=0; e<N_e; e++)
        file << "4\t" << Elements[e][0] << "\t" << Elements[e][1] << "\t" << Elements[e][2] << "\t" << Elements[e][3] << endl;
    
    file << "CELL_TYPES " << N_e << endl;
    for(int e=0; e<N_e; e++)
        file << "10" << endl;
	
    file << "POINT_DATA " << N_p << endl;
    file << "SCALARS phi double 1" << endl;
    file << "LOOKUP_TABLE \"default\"" << endl;
    for(int p=0; p<N_p; p++)
        file << setprecision(5) << T[p] << endl;
    file.close();
	
	return;
}














