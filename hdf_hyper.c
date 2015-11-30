// From x-y ascii data from xgrafix, reconstruct the density profile
#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "string.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

#define field	0
#define density	1
#define raman	2



main(int argc, char *argv[])
{
   char fileName[100],dataName[100],name[100];
   int dimension,step,M,N,xinter,yinter,zinter;
   int i,j,k,nx,ny,nz,fileType,offset[3];
   int minX,minY,minZ,maxX,maxY,maxZ,remainY,remainZ,tmp,subY,subZ;
   int nySub,nzSub,rank,rankY,rankZ,minYSub,maxYSub,minZSub,maxZSub;
   int startI,startJ,startK,startY,startZ,centerI,centerJ,centerK;
   int nxPrime,nyPrime,nySubPrime,nzSubPrime;
   float ***data;
   void restoreIntMeta();
   int whatDataType();
   void restoreData2D();

   int myrank, nTasks;
   MPI_Status status;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(argc < 9)
   {
      printf("hdf_hyper dimension step fileType dataName L M xinterval yinterval zinterval\n");
      exit(0);
   }

   dimension=atoi(argv[1]);
   step=atoi(argv[2]);
   fileType=whatFileType(argv[3]);
//   dataName=argv[4];
   M=atoi(argv[5]);
   N=atoi(argv[6]);
   xinter=atoi(argv[7]);
   yinter=atoi(argv[8]);
   zinter=atoi(argv[9]);

   switch (fileType)   {
   case field :
     sprintf(fileName,"field%d.h5",step);
     break;
   case density :
     sprintf(fileName,"density%d.h5",step);
     break;
   case raman :
     sprintf(fileName,"raman%d.h5",step);
     break;
   }

   if(dimension>1)
   {
     sprintf(dataName,"/nx");   
     restoreIntMeta(fileName,dataName,&nx);
     sprintf(dataName,"/ny");   
     restoreIntMeta(fileName,dataName,&ny);
     nz=1;
   }
   if(dimension>2)
   {
     sprintf(dataName,"/nz");   
     restoreIntMeta(fileName,dataName,&nz);
   }

   //assigning to ranks
   nySub=ny/M;
   subY=nySub;
   remainY=ny%M;
   nzSub=nz/N;
   subZ=nzSub;
   remainZ=nz%N;
   for(rankZ=0; rankZ<N; rankZ++)
   {
     minY=maxY=0;
     for(rankY=0; rankY<M; rankY++)
     {
       rank=rankY+(rankZ*M);
       if(rankY<remainY)   tmp=subY+1;
       else                tmp=subY;
       minY=maxY;
       maxY=minY+tmp;
       if(myrank==rank)
       {
         nySub=tmp;
         minYSub=minY;
         maxYSub=maxY;
       }
     }  
   }
   for(rankY=0; rankY<M; rankY++)
   {
     minZ=maxZ=0;
     for(rankZ=0; rankZ<N; rankZ++)
     {
       rank=rankY+(rankZ*M);
       if(rankZ<remainZ)   tmp=subZ+1;
       else                tmp=subZ;
       minZ=maxZ;
       maxZ=minZ+tmp;
       if(myrank==rank)
       {
          minZSub=minZ;
          maxZSub=maxZ;
          nzSub=tmp;
       }
     }
   }
  
   data=(float ***)malloc(nx*sizeof(float **));
   for(i=0; i<nx; i++)
   {
     data[i]=(float **)malloc(nySub*sizeof(float *));
     for(j=0; j<nySub; j++)
       data[i][j]=(float *)malloc(nzSub*sizeof(float ));
   }
   
   switch (fileType)   {
   case field :
     if(dimension==2)
     {
       offset[0]=0;
       offset[1]=minYSub;
//here     
       restoreData2D(data,fileName,argv[4],nx,ny,nySub,0,nySub,offset);
       centerI=(int)(0.5*nx);
       startI=centerI%xinter;
       nxPrime=(int)((nx-startI+xinter-1)/xinter);
       centerJ=(int)(0.5*ny);
       startJ=centerJ%yinter;
       nyPrime=(int)((ny-startJ+yinter-1)/yinter);
       tmp=0;
       while((minYSub-startJ+tmp)%yinter!=0)
         tmp++;
       startY=minYSub+tmp;
       nyPrime=(int)((maxYSub-startY+yinter-1)/yinter);
printf("myrank=%d, centerJ=%d, startJ=%d,startY=%d,nyPrime=%d\n",myrank,centerJ,startJ,startY,nyPrime);
     }

     break;
   case density :
     if(dimension==2)
     {
       offset[0]=0;
       offset[1]=minYSub;
     
       sprintf(dataName,"/0");
       restoreData2D(data,fileName,dataName,nx,ny,nySub,0,nySub,offset);


     }    
     break;
   }
      
//   else if(fileName!=NULL && dataName==NULL)
//     nx=atoi(argv[3]);
//   else
//   {
//     printf("file name is missing!\n");
//     exit(0);
//   }  


   MPI_Finalize();


//      rangeY=atof(argv[7]);
//      ymax = ymin+rangeY;

/*
   ne=(float **)malloc((Nx+1)*sizeof(float *));
   xc=(float *)malloc((Nx+1)*sizeof(float));
   yc=(float *)malloc((Ny+1)*sizeof(float));
   for(i=0;i<=Nx;i++)  ne[i]=(float *)malloc((Ny+1)*sizeof(float));

   for(i=0;i<=Nx;i++)
      for(j=0;j<=Ny;j++)  ne[i][j]=0.0;

   for(i=0;i<=Nx;i++)  xc[i]=xmin+i*dx;
   for(j=0;j<=Ny;j++)  yc[j]=ymin+j*dy;

   if(mode==1 || mode==0)
   {
     dataX=(float *)malloc((dataNum)*sizeof(float));
     dataY=(float *)malloc((dataNum)*sizeof(float));
     sprintf(fileName,"%dParticle%d.h5",species1,step);

     sprintf(name,"/x");
     restoreFloatArray(fileName,name,dataX,dataNum);
     sprintf(name,"/y");
     restoreFloatArray(fileName,name,dataY,dataNum);
     for(ii=0; ii<dataNum; ii++)
     {
       i=(int)((dataX[ii]-xmin)/dx);
       j=(int)((dataY[ii]-ymin)/dy);
       if(0 < i && i< Nx && 0 < j && j < Ny)
       {
         ne[i][j]+=(xc[i+1]-dataX[ii])*(yc[j+1]-dataY[ii])/dx/dy*np2c1;
         ne[i+1][j]+=(dataX[ii]-xc[i])*(yc[j+1]-dataY[ii])/dx/dy*np2c1;
         ne[i][j+1]+=(xc[i+1]-dataX[ii])*(dataY[ii]-yc[j])/dx/dy*np2c1;
         ne[i+1][j+1]+=(dataX[ii]-xc[i])*(dataY[ii]-yc[j])/dx/dy*np2c1;
       }
     }
     free(dataX);
     free(dataY);
*/
/*
   for (i=0; i<=Nx;i++)
   {
      for(j=0; j<=Ny; j++)
      {
         printf("%g  %g  %g\n",xc[i],yc[j],ne[i][j]/dx/dy);
      }
      printf("\n");
   }
*/
}

void restoreData2D(float ***data,char *fileName,char *dataName,int nx,int ny,int nySub,int jstart,int jend,int *offSet)
{
  int i,j,k,start;
  float *memory;
  char name[100];
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[3],count[3],offset[3];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=nx;
  dimsf[1]=ny;
  filespace=H5Screate_simple(2,dimsf,NULL);

  count[1]=nySub;
  count[0]=nx;
  offset[1]=offSet[1];
  offset[0]=offSet[0];
  memspace=H5Screate_simple(2,count,NULL);

  memory = (float *)malloc(nx*nySub*sizeof(float ));
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_NATIVE_FLOAT,memspace,subfilespace,plist_id,memory);
/*
  start=0;
  for(i=0; i<nx; i++)
  {
    for(j=jstart; j<jend; j++)
      data[i][j][0]=memory[start+j-jstart];
    start+=nySub;
  }
*/
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
  free(field);
}












void restoreIntMeta(char *fileName,char *dataName,int *data)
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
        
void restoreFloatArray(char *fileName,char *dataName,float *data,int totalCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=totalCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

int whatFileType(char *str)
{
   if(strstr(str,"field"))           return field;
   else if(strstr(str,"density"))      return density;
   else if(strstr(str,"raman"))    return raman;
   else return 0;
}

