#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"
 
void main(int argc, char *argv[])
{
  float minX,maxX,minE,maxE,dx,eenergy,superP;
  int i,j,rnk,numParticle,numX,numE,time,indexI,indexJ,index;
  int cnt,min,max,remain,tmp,numData;
  float **n,*nE,*x,*E,*energy,*xx,*px,*py,*pz,*send,*sendE;
  FILE *in,*out;
  char fileName[100];
  void restoreData();

  MPI_Init(&argc,&argv);

  int myrank, nTasks;    
  MPI_Status status;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  int cntOffSet[nTasks];
  for(i=0; i<nTasks; i++)
    cntOffSet[i]=0;

  if(argc < 8)   { 
    printf("energy3D_HDF numParticle minX maxX numX minE maxE time superP(3.14e4)\n");
    exit(0);  
  }

  numParticle=atoi(argv[1]);
  minX=atof(argv[2]);
  maxX=atof(argv[3]);
  numX=atoi(argv[4]);
  minE=atof(argv[5]);
  maxE=atof(argv[6]);
  time=atoi(argv[7]);
  superP=atof(argv[8]);
  

  dx=(maxX-minX)/((float) numX);
  numE=(int)((maxE-minE)/1e6);
  x=(float *)malloc(numX*sizeof(float ));
  E=(float *)malloc(numE*sizeof(float ));
  nE=(float *)malloc(numE*sizeof(float ));
  n=(float **)malloc(numX*sizeof(float *));
  for(i=0; i<numX; i++)
    n[i]=(float *)malloc(numE*sizeof(float ));
  for(i=0; i<numX; i++)
    for(j=0; j<numE; j++)
      n[i][j]=0.0;
  for(i=0; i<numE; i++)
    nE[i]=0.0;

  cnt=numParticle/nTasks;
  remain=numParticle%nTasks;
  min=max=0;
  for(i=0; i<nTasks; i++)
  {
    if(i<remain) tmp=cnt+1;
    else	 tmp=cnt;
    min=max;
    max=min+tmp;
    if(i==myrank)
    {
      cntOffSet[i]=min;
      numData=tmp;
    }
  }
  xx=(float *)malloc(numData*sizeof(float ));
  px=(float *)malloc(numData*sizeof(float ));
  py=(float *)malloc(numData*sizeof(float ));
  pz=(float *)malloc(numData*sizeof(float ));
  energy=(float *)malloc(numData*sizeof(float ));
  send=(float *)malloc(numX*numE*sizeof(float ));
  sendE=(float *)malloc(numE*sizeof(float ));

  sprintf(fileName,"0Particle%d.h5",time);
  restoreData(xx,fileName,"/x",numParticle,numData,cntOffSet);
  restoreData(px,fileName,"/px",numParticle,numData,cntOffSet);
  restoreData(py,fileName,"/py",numParticle,numData,cntOffSet);
  restoreData(pz,fileName,"/pz",numParticle,numData,cntOffSet);

  //Saving denE file
  for(i=0; i<numData; i++)
  {
    indexI=(int)((xx[i]-minX)/dx);
    eenergy=(sqrt(1.0+px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])-1.0)*0.511*1e6;
    indexJ=(int)((eenergy-minE)/1e6);
    if(indexI>=0 && indexI<numX && indexJ>=0 && indexJ<numE)
      n[indexI][indexJ]+=1.0;
  }
  index=0;
  for(i=0; i<numX; i++)
    for(j=0; j<numE; j++)
    {
      send[index]=n[i][j];
      index++;
    }
  
  for(rnk=1; rnk<nTasks; rnk++)
  {
    if(myrank==rnk)
      MPI_Send(send,numE*numX,MPI_FLOAT,0,myrank,MPI_COMM_WORLD);
    else if(myrank==0)
    {
      MPI_Recv(send,numE*numX,MPI_FLOAT,rnk,rnk,MPI_COMM_WORLD,&status);
      index=0;
      for(i=0; i<numX; i++)
        for(j=0; j<numE; j++)
        {
          n[i][j]+=send[index];
          index++;
        }
    }
    else	;
    printf("rank%d is done\n",rnk);
  }

  for(i=0; i<numX; i++)
    x[i]=i*dx+minX;
  for(i=0; i<numE; i++)
    E[i]=i*1e6+minE;

  if(myrank==0)
  {
    sprintf(fileName,"denE%d",time);
    out=fopen(fileName,"w");
    for(i=0; i<numX; i++)
    {
      for(j=0; j<numE; j++)
        fprintf(out,"%g %g %g\n",x[i],E[j],n[i][j]*superP*1.6e-19*1e12);
      fprintf(out,"\n");
    }
    fclose(out);
    printf("%s is made.\n",fileName);
  }  
  else	;

  //Saving energy file
  for(i=0; i<numData; i++)
  {
    indexI=(int)((xx[i]-minX)/dx);
    eenergy=(sqrt(1.0+px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])-1.0)*0.511*1e6;
    indexJ=(int)((eenergy-minE)/1e6);
    if(indexJ>=0 && indexJ<numE)
      nE[indexJ]+=1.0;
  }
  index=0;
  for(j=0; j<numE; j++)
  {
    sendE[index]=nE[j];
    index++;
  }
  
  for(rnk=1; rnk<nTasks; rnk++)
  {
    if(myrank==rnk)
      MPI_Send(sendE,numE,MPI_FLOAT,0,myrank,MPI_COMM_WORLD);
    else if(myrank==0)
    {
      MPI_Recv(sendE,numE,MPI_FLOAT,rnk,rnk,MPI_COMM_WORLD,&status);
      index=0;
      for(j=0; j<numE; j++)
      {
        nE[j]+=sendE[index];
        index++;
      }
    }
    else	;
    printf("rank%d is done\n",rnk);
  }

  for(i=0; i<numE; i++)
    E[i]=i*1e6+minE;

  if(myrank==0)
  {
    sprintf(fileName,"energy%d",time);
    out=fopen(fileName,"w");
    for(j=0; j<numE; j++)
      fprintf(out,"%g %g\n",E[j],nE[j]*superP*1.6e-19*1e12);
    fclose(out);
    printf("%s is made.\n",fileName);
  }  
  else	;

  free(x);
  free(E);
  free(xx);
  free(px);
  free(py);
  free(pz);
  free(energy);
  free(send);
  free(sendE);
  for(i=0; i<numX; i++)
    free(n[i]);
  free(n);
  free(nE);

  MPI_Finalize();
}


void restoreData(float *data,char *fileName,char *dataName,int totalCnt,int cnt,int *cntOffSet)
{
  hid_t file_id,dset_id,plist_id,group_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[1],count[1],offset[1];

  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//  H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//  MPI_Barrier(MPI_COMM_WORLD);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=totalCnt;     
  filespace=H5Screate_simple(1,dimsf,NULL);

  count[0]=cnt;
  offset[0]=cntOffSet[myrank];
  memspace=H5Screate_simple(1,count,NULL);

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void saveMetaData(char *fileName,char *dataName,int *data)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

