//////////////////////////////////////////////////////////////////////////
//
// Cluster and Cloud Computing Assignment 1
//
// Student Name: Jie JIN
//
// Student No:   652600
//
// Problem:		 HPC Data Processing
//
// Compilation:  mpicxx -o run  ass1.cpp 
//
// Execution:	 mpirun -n NUM_PROCESSORS ./run -f FILE_NAME -k KEYWORD
//
// Edward:       qsub run.pbs
//				
/////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <mpi.h>
using namespace std;
#define NOTFOUND ""
#define TAG 1
#define LBOUND 0
#define RBOUND 1
#define ROOT 0
#define ERROR -1
#define P(VAR) cout << #VAR " = " << VAR << endl



void ErrorMessage(int error, int worldRank, string s) {
	cerr<< "Process " << worldRank;
	error == ERROR ? cerr << " " : cerr << " Error " << error;
	cerr << " in " << s  <<endl;
	MPI_Finalize();
	exit(-1);
}

string GetArg(int argc, char *argv[], string arg) {
    for (int i = 1; i < argc - 1; ++i)
        if (strcmp(argv[i], arg.c_str()) == 0)
			return (argv[i+1]); 
	return NOTFOUND;
}

MPI_Offset GetSize(string fileName, int worldRank) {
	int error;
	MPI_Offset  fileSize;
	MPI_File file;
	error = MPI_File_open(MPI_COMM_WORLD, (char*)fileName.c_str(),
                  MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	if(error != MPI_SUCCESS)
		ErrorMessage(error, worldRank, "MPI_File_open()");
	error = MPI_File_get_size(file, &fileSize);
	if(error != MPI_SUCCESS)
		ErrorMessage(error, worldRank, "MPI_File_get_size()");
	MPI_File_close(&file);
	return fileSize;
}

vector<MPI_Offset> EstimateOffsets(MPI_Offset size,
								   const int numOfNodes) {
	vector<MPI_Offset> offsets(numOfNodes << 1, 0);
	if (numOfNodes != 1) {
		MPI_Offset quotient = size / numOfNodes;
		for (size_t i = 0; i != (size_t)numOfNodes - 1; 
			 offsets[i*2+1] = (i+1)* quotient, ++i);
	}
	offsets.back() = size - 1; 

	return offsets;
}

void AdjustOffsets(vector<MPI_Offset> &offsets,
				   string fileName, int worldRank) {
	char ch;
	ifstream file(fileName.c_str());
	
	if (file.is_open()) {
		for (size_t i = 1; i !=offsets.size() -1; i += 2) {
			off_t tmp = offsets[i];
			while (true){
				file.seekg(tmp, file.beg);
				if ((file >> noskipws >> ch) && ch != ' '){
					tmp++;
				} else {
					offsets[i] = tmp-1;
					offsets[i+1] = offsets[i] + 2;
					break;
				}
			}
		}	
		file.close();
	} else {
		ErrorMessage(ERROR, worldRank, "ifstream file()");
	}

}

int isNotAlphaNum(char c){
	return !isalnum(c);
}

void StripCreateTmpFile(string fileName, int worldRank,
						vector<MPI_Offset> range) {
	ifstream file(fileName.c_str());
	char tmpFileName[20];
	snprintf(tmpFileName, 10, "%d", worldRank);
	ofstream tmpFile(tmpFileName);
	string b;
	
	if (file.is_open()){
		off_t length = range[RBOUND] - range[LBOUND] + 1;
		file.seekg(range[LBOUND]);
		char *buffer = new char[length];
		file.read(buffer, length);
		b = string(buffer);
		replace_if(b.begin(), b.end(), isNotAlphaNum, ' ');
		transform(b.begin(), b.end(), b.begin(), ::tolower);
		tmpFile.write(b.c_str(), length);
		file.close();
		tmpFile.close();
		delete[] buffer;
	} else {
		ErrorMessage(ERROR, worldRank, "ifstream file()");
	}
}

long ScanTmpFile(int worldRank, string keyword) {
	char fileName[20];
	snprintf(fileName, 20, "%d", worldRank);
	ifstream file(fileName);
	long result = 0;
	string token;
	
	if (file.is_open()) {
		while(!file.eof())
			while(file >> token)
				if (token == keyword)
					++result;
	} else {
		ErrorMessage(ERROR, worldRank, "ifstream file()");
	}
	
	return result;
}

void DeleteTmpFile(int worldRank) {
	char fileName[20];
	snprintf(fileName, 20, "%d", worldRank);
	if(remove(fileName) != 0)
		ErrorMessage(ERROR, worldRank, "remove()");
}


void OutputResult(double wtime, long result,
				  string keyword, int worldSize){
	cout << keyword <<" #" << result << endl;
	wtime	= MPI_Wtime() - wtime;
	cout << "Simulation took " << wtime << " seconds with "
		 << worldSize << " processes" << endl;
}

int main(int argc, char *argv[]) {
	int worldRank, worldSize, error;
	double wtime = 0;
	long  result = 0, tmp = 0;
	MPI_Offset fileSize;
	string fileName, keyword;
	vector<MPI_Offset>  offsets;
	vector<MPI_Offset>  range(2, 0);
	MPI_Status status;
	MPI_Datatype rangeT;
	
	//Initialize MPI environment
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
	
	//create new MPI data type for partition the range of each file
	MPI_Type_contiguous(2, MPI_LONG_LONG, &rangeT);
	MPI_Type_commit(&rangeT);
	
	//Get command line
	if (argc != 5) {
		if (worldRank == ROOT)
			cerr << "Usage: mpirun -np <numProcs> "
				 << argv[0] << " -f <fileName> -k <keyword>" << endl;
		MPI_Finalize();
		exit(-1);
	}
	fileName = GetArg(argc, argv, "-f");
	keyword = GetArg(argc, argv, "-k");

	//Strip keyword
	int kCount = count_if (keyword.begin(), keyword.end(), isNotAlphaNum);
	if (kCount != 0) {
		if (worldRank == ROOT)
			cerr << "Error: Keyword should not "
				 << "contain any hyphenated input data."<< endl;
		MPI_Finalize();
		exit(-1);
	}
	transform(keyword.begin(), keyword.end(), keyword.begin(), ::tolower);

	//Each processor get file size
	fileSize= GetSize(fileName,worldRank);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (worldRank == ROOT){
        // ROOT processor calculate offsets and partition them
		wtime	= MPI_Wtime();
		offsets = EstimateOffsets(fileSize, worldSize);
		if (worldSize > 1) {
			AdjustOffsets(offsets, fileName, worldRank);
			for (int i = 1 , j = 2; i != worldSize;
				 error = MPI_Send(&offsets[j], 1, rangeT, i++,
								  TAG, MPI_COMM_WORLD), j+=2);
			if(error != MPI_SUCCESS)
				ErrorMessage(error, worldRank, "MPI_Send()");
		}
		range[LBOUND] = offsets[0]; range[RBOUND] = offsets[1];
	} else {
		if (worldSize > 1) {
			error = MPI_Recv(&range[0], 1, rangeT, ROOT, TAG,
							 MPI_COMM_WORLD, &status);
			if(error != MPI_SUCCESS)
				ErrorMessage(error, worldRank, "MPI_Recv()");
		}
	}

	// Strip the origional file and create  smaller temporary files
	StripCreateTmpFile(fileName, worldRank, range);
	// Count the number of times a given term (word/string) appears. 
	result = ScanTmpFile(worldRank, keyword);
	// Delete temp files
	DeleteTmpFile(worldRank);
	
	// send results to ROOT node.
	// ROOT node outputs result
	MPI_Barrier(MPI_COMM_WORLD);
	if (worldSize > 1){
		if (worldRank == ROOT){
			for (int i = 1; i != worldSize; result += tmp) {
				error = MPI_Recv(&tmp, 1, MPI_LONG, i++, MPI_ANY_TAG,
								 MPI_COMM_WORLD,&status);
				if(error != MPI_SUCCESS)
					ErrorMessage(error, worldRank, "MPI_Recv()");
			}
			OutputResult(wtime, result, keyword, worldSize);
		} else {
			error = MPI_Send(&result, 1, MPI_LONG, ROOT, worldRank,
							 MPI_COMM_WORLD);
			if(error != MPI_SUCCESS)
				ErrorMessage(error, worldRank, "MPI_Send()");
		}
	} else {
		OutputResult(wtime, result, keyword, worldSize);
	}

	MPI_Finalize();
}
