#include "utility.h"

int main (int argc, char **argv) {

	MPI_Init(&argc, &argv);
	int world_rank;
	int world_size;
	int I, J;
	MPI_Status status;
	MPI_Datatype row, col, colb, all;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int Pi = getIntArg(argc, argv, "-i", 2);
	int Pj = getIntArg(argc, argv, "-j", 2);
	I = MAPI / Pi;
	J = MAPJ / Pj;
	int ****map =  alloc4d(world_size,I + 2, J + 2 , 4);
	int ****col_buffer;
	int ****row_buffer;
	int neighbours[4] = {getL(Pj), getR(Pj), getT(Pj), getB(Pj)};
	double wtime;
	
	MPI_Type_vector(J , 4, 4   , MPI_INT, &row);
	MPI_Type_commit(&row);
	MPI_Type_vector(I , 4, 4*(J+2) , MPI_INT, &col);
	MPI_Type_commit(&col);
	MPI_Type_vector(I , 4, 4 , MPI_INT, &colb);
	MPI_Type_commit(&colb);
	MPI_Type_vector(I , 1, 2*sizeof(row)/J + 1, row, &all);
	MPI_Type_commit(&all);
	if(world_rank==0)
		wtime	= MPI_Wtime();
	initMap(map[world_rank], I, J);

	row_buffer = alloc4d(1,1, J ,4);
	col_buffer = alloc4d(1,I, 1, 4);

	for (int k = 0; k <TIME; k++ ) {
		//send the top side to the upper processor
		if ((world_rank / Pj) % 2 == 0) {
			MPI_Send(&map[world_rank][1][1][0], 1, row, neighbours[UP], 1, MPI_COMM_WORLD);
		} else{
			MPI_Recv(&map[world_rank][I+1][1][0] , 1, row, neighbours[DOWN],1, MPI_COMM_WORLD, &status );
			MPI_Send(&map[world_rank][1][1][0], 1, row, neighbours[UP], 2, MPI_COMM_WORLD);}
		if ((world_rank / Pj) % 2 == 0) {
			MPI_Recv( &map[world_rank][I+1][1][0] , 1, row, neighbours[DOWN], 2,MPI_COMM_WORLD, &status );
		}

		//send the bottom side to the bottom processor
		if ((world_rank / Pj) % 2 == 0) {
			MPI_Send(&map[world_rank][I][1][0], 1, row, neighbours[DOWN], 0,MPI_COMM_WORLD);
		} else{
			MPI_Recv(&map[world_rank][0][1][0] , 1, row, neighbours[UP], 0, MPI_COMM_WORLD, &status);
			MPI_Send(&map[world_rank][I][1][0], 1, row, neighbours[DOWN], 0,MPI_COMM_WORLD);
		}
		if ((world_rank / Pj) % 2 == 0) {
			MPI_Recv(&map[world_rank][0][1][0] , 1, row, neighbours[UP], 0,  MPI_COMM_WORLD, &status);
		}

		//send the left side to the left processor
		if (world_rank % 2 == 0) {
			MPI_Send(&map[world_rank][1][1][0], 1, col, neighbours[LEFT], 5, MPI_COMM_WORLD);
		} else {
			MPI_Recv(&map[world_rank][1][J+1][0] , 1, col, neighbours[RIGHT], 5, MPI_COMM_WORLD, &status);
			MPI_Send(&map[world_rank][1][1][0], 1, col, neighbours[LEFT], 6, MPI_COMM_WORLD);
		}
		if (world_rank %2 == 0) {
			MPI_Recv(&map[world_rank][1][J+1][0] , 1, col, neighbours[RIGHT], 6, MPI_COMM_WORLD, &status);
		}
		//send the right side to the right processor 
		if (world_rank % 2 == 0) {
			MPI_Send(&map[world_rank][1][J][0], 1, col, neighbours[RIGHT], 7, MPI_COMM_WORLD);
		} else {
			MPI_Recv(&map[world_rank][1][0][0] , 1, col, neighbours[LEFT], 7, MPI_COMM_WORLD, &status);
			MPI_Send(&map[world_rank][1][J][0], 1, col, neighbours[RIGHT], 8, MPI_COMM_WORLD);
		}
		if (world_rank %2 == 0) {
			MPI_Recv(&map[world_rank][1][0][0] , 1, col, neighbours[LEFT], 8, MPI_COMM_WORLD, &status);
		}

	   
		MPI_Barrier(MPI_COMM_WORLD);
		move(map[world_rank],I,J);
		MPI_Barrier(MPI_COMM_WORLD);

		// send the upper boundary back to its top processor
		if ((world_rank / Pj) % 2 == 0) {
			MPI_Send(&map[world_rank][0][1][0], 1, row, neighbours[UP], 9,MPI_COMM_WORLD);
		} else {
			MPI_Recv(&row_buffer[0][0][0][0] , 1, row, neighbours[DOWN], 9, MPI_COMM_WORLD, &status);
			MPI_Send(&map[world_rank][0][1][0], 1, row, neighbours[UP], 10,MPI_COMM_WORLD);
		}
		if ((world_rank / Pj) % 2 == 0) {
	  	  MPI_Recv(&row_buffer[0][0][0][0] , 1, row, neighbours[DOWN], 10, MPI_COMM_WORLD, &status);
		}
	
		decisionForBoundary(row_buffer[0], map[world_rank], J, I, I, J);

		// send the bottom boundary back to its bottom processor
		if ((world_rank / Pj) % 2 == 0) {
			MPI_Send(&map[world_rank][I+1][1][0], 1, row, neighbours[DOWN], 11,MPI_COMM_WORLD);
		} else {
		    MPI_Recv(&row_buffer[0][0][0][0] , 1, row, neighbours[UP], 11, MPI_COMM_WORLD, &status);
			MPI_Send(&map[world_rank][I+1][1][0], 1, row, neighbours[DOWN], 12,MPI_COMM_WORLD);
		}
		if ((world_rank / Pj) % 2 == 0) {
			MPI_Recv(&row_buffer[0][0][0][0] , 1, row, neighbours[UP], 12, MPI_COMM_WORLD, &status);
		}
		decisionForBoundary(row_buffer[0], map[world_rank], J, 1, I, J);
	
		// send the left boundary back to its left processor
		if (world_rank % 2 == 0) {
			MPI_Send(&map[world_rank][1][0][0], 1, col, neighbours[LEFT], 13, MPI_COMM_WORLD);
		} else {
			MPI_Recv(&col_buffer[0][0][0][0] , 1, colb, neighbours[RIGHT], 13, MPI_COMM_WORLD, &status);
			MPI_Send(&map[world_rank][1][0][0], 1, col, neighbours[LEFT], 14, MPI_COMM_WORLD);
		}
		if (world_rank %2 == 0) {
			MPI_Recv(&col_buffer[0][0][0][0] , 1, colb, neighbours[RIGHT], 14, MPI_COMM_WORLD, &status);;
		}
		decisionForBoundary(col_buffer[0], map[world_rank], I, J, I, J);

	
		// send the right boundary back to its right processor
		if (world_rank % 2 == 0) {
			MPI_Send(&map[world_rank][1][J+1][0], 1, col, neighbours[RIGHT], 15, MPI_COMM_WORLD);
		} else {
			MPI_Recv(&col_buffer[0][0][0][0] , 1, colb, neighbours[LEFT], 15, MPI_COMM_WORLD, &status);
			MPI_Send(&map[world_rank][1][J+1][0], 1, col, neighbours[RIGHT], 16, MPI_COMM_WORLD);
		}
		if (world_rank %2 == 0) {
			MPI_Recv(&col_buffer[0][0][0][0] , 1, colb, neighbours[LEFT], 16, MPI_COMM_WORLD, &status);
		}
		decisionForBoundary(col_buffer[0], map[world_rank], I, 1, I, J);
		//update map & reset property
		updateMap(map[world_rank], I, J);
	  
		if (world_rank!=0) {
			for (int i = 1; i< world_size; i++){
				if (world_rank == i)
					MPI_Send(&map[world_rank][1][1][0],   1, all, 0, world_rank, MPI_COMM_WORLD);
			}
		}
		if (world_rank==0)  {
			for (int i = 1; i< world_size; i++){
				MPI_Recv(&map[i][1][1][0],   1, all, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}		  
		} 

		if (world_rank == 0){
				calPup(map, I, J);
		}
	}

	if(world_rank == 0) {
		wtime	= MPI_Wtime() - wtime;
		cout << "Simulation took " << wtime << " seconds with "
			 << world_size << " processes" << endl;
	}


	MPI_Finalize();

}

