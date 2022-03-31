#include "mpi.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>
#include "random.h"

using namespace std;



int main(int argc, char* argv[]){

  int size, rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status stat;
  vector<int> imesg1(6),imesg2(6), idx={0,1};
  int N=6;
//  if(rank==0) imesg1 = {0,0,0,0,0,0};
//  if(rank==1) imesg2 = {1,2,3,4,5,6};

  for(int j=0; j<1; j++){
    int tag1 = j*2;
    int tag2 = j*2+1;
    if(rank == idx[2*j]){
      imesg1.clear();
      imesg1 = {6,5,4,3,2,1};
      MPI_Send(&imesg1.front(), N, MPI_INT, idx[2*j+1], tag1, MPI_COMM_WORLD);
      MPI_Recv(&imesg2.front(), N, MPI_INT, idx[2*j+1], tag2, MPI_COMM_WORLD, &stat);
      for(unsigned int i=0; i<6; i++) cout << imesg2[i] << " ";
      cout << endl;
    }
    if(rank == idx[2*j+1]){
      imesg2.clear();
      imesg2 = {1,2,3,4,5,6};
      MPI_Recv(&imesg1.front(), N, MPI_INT, idx[2*j], tag1, MPI_COMM_WORLD, &stat);
      MPI_Send(&imesg2.front(), N, MPI_INT, idx[2*j], tag2, MPI_COMM_WORLD);
      for(unsigned int i=0; i<6; i++) cout << imesg1[i] << " ";
      cout << endl;
    }
  }
/*
  if(rank==0){
    for(unsigned int i=0; i<6; i++) cout << imesg1[i] << " ";
    cout << endl;
  }
*/
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();


return 0;
}


//if(rank==1) MPI_Send(&imesg.front(),6,MPI_INT,0,itag,MPI_COMM_WORLD);
//else if(rank==0)  MPI_Recv(&imesg.front(),6,MPI_INT,1,itag,MPI_COMM_WORLD,&stat);
//MPI_Bcast(my_values,3,MPI_INTEGER,0, MPI_COMM_WORLD);
//MPI_Gather(&isend,1,MPI_INTEGER,irecv,1,MPI_INTEGER,0,MPI_COMM_WORLD);
