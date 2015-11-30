#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hdf5.h"
#include "hdf5_hl.h"
 

void main(int argc, char *argv[])
{
   int i,j,k,nx,ny,nz,time;
   float *dataX,*dataY,*dataZ;
   float dt,kp,v0,gamma0;
   char fileName[100],dataName[100];
   void readCoord();
   void writeCoord();

   MPI_Init(&argc,&argv);

   int myrank, nTasks;    
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   hid_t file_id,group_id,dset_id,filespace;
   hsize_t metaDim[1],dimsf[3],coord[1];
   herr_t status;

   if(argc<4)  {
      printf("hdfConvert time nx ny nz\n");
      exit(0);
   }

   time=atoi(argv[1]);
   nx=atoi(argv[2]);
   ny=atoi(argv[3]);
   nz=atoi(argv[4]);
 
   kp=4.21e5;
   dt=1.67e-16;
   gamma0=10;
   v0=3e8*(1.0-0.5/gamma0/gamma0);

   sprintf(fileName,"field%d.h5",time);
   dataX=(float *)malloc(nx*sizeof(float ));
   dataY=(float *)malloc(ny*sizeof(float ));
   dataZ=(float *)malloc(nz*sizeof(float ));
   readCoord(dataX,fileName,"X",nx,time);
   readCoord(dataY,fileName,"Y",ny,time);
   readCoord(dataZ,fileName,"Z",nz,time);
   for(i=0; i<nx; i++)
     dataX[i]=dataX[i]*kp-time*dt*v0*kp;
   for(i=0; i<ny; i++)
     dataY[i]=dataY[i]*kp;
   for(i=0; i<nz; i++)
     dataZ[i]=dataZ[i]*kp;
   writeCoord(dataX,fileName,"X",nx,time);
   writeCoord(dataY,fileName,"Y",ny,time);
   writeCoord(dataZ,fileName,"Z",nz,time);
//   for(i=0; i<nx; i++)
//     printf("datsX[%d]=%g\n",i,dataX[i]);
/*
    metaDim[0]=1;
    file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);
    filespace=H5Screate_simple(1,metaDim,NULL);
    dset_id=H5Dopen2(file_id,"/minXSub",H5P_DEFAULT);
    status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,&(D->minXSub));
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Fclose(file_id);

    else	;
    MPI_Bcast(&(D->minXSub),1,MPI_INT,0,MPI_COMM_WORLD);
printf("myrank=%d, minXSub=%d\n",myrank,D->minXSub);
  
    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Split-1)*3+2:
      ny=D->ny+5;      
      nx=D->nx+5;     
      istart=0;
      iend+=3; 
      jstart=0;
      jend+=3;
      nxSub+=5;
      nySub+=5;
      offSetY=D->minYSub-D->minYDomain;
      
//      MPI_Gather(&offSetY,1,MPI_INT,offSetJ,1,MPI_INT,0,MPI_COMM_WORLD);
      restoreFileField2D(D->Ex,"/Ex",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->Pr,"/Pr",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->Pl,"/Pl",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->Bx,"/Bx",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->Sr,"/Sr",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->Sl,"/Sl",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->ExC,"/ExC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->PrC,"/PrC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->PlC,"/PlC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->BxC,"/BxC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->SrC,"/SrC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->SlC,"/SlC",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->Jx,"/Jx",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->Jy,"/Jy",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->Jz,"/Jz",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->JxOld,"/JxOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->JyOld,"/JyOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
      restoreFileField2D(D->JzOld,"/JzOld",nx,ny,nxSub,nySub,istart,iend,jstart,jend,iteration,D->dimension,offSetY);
*/

   MPI_Finalize();
}


void readCoord(float *data,char *fileName,char *dataName,int numData,int time)
{
  int i,j,k,start;
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,group_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t metaDim[1];

  metaDim[0]=numData;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void writeCoord(float *data,char *fileName,char *dataName,int numData,int time)
{
  int i,j,k,start;
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,group_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t metaDim[1];

  metaDim[0]=numData;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}
