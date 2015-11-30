// From x-y ascii data from xgrafix, reconstruct the density profile
#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "math.h"

void restoreFloatArray(char *fileName,char *dataName,float *data,int totalCnt);
void saveIntMeta(char *fileName,char *dataName,int *data);
void saveDataHDF3D(char *fileName,float ***ne,int numX,int numY,int numZ,float *xc,float *yc,float *zc);
void saveXmf(char *fileName,int numX,int numY,int numZ);
void restoreHyperSlab2D(char *fileName,char *dataName,float *data,int nx,int ny,int pick);
void restoreField2D(char *fileName,char *dataName,float **data,int nx,int ny);

main(int argc, char *argv[])
{
   FILE *in,*out;
   char fileName[100],dataName[100];
   int i,j,nx,ny,pick,mode,dimension,stX,stY;
   float pickY;
   float **Ex,*dataX,*dataY;

   if(argc < 4)
   {
      printf("hdf_field mode dimension \n");
      printf("mode:0 => fileName(.h5) dataName strideX strideY\n");
      exit(0);
   }
   mode=atoi(argv[1]);
   dimension=atoi(argv[2]);

   if(mode==0 && dimension==2) 
   {
      sprintf(fileName,"%s.h5",argv[3]);
      sprintf(dataName,"%s",argv[4]);
      stX=atoi(argv[5]);
      stY=atoi(argv[6]);

      saveIntMeta(fileName,"nx",&nx);
      saveIntMeta(fileName,"ny",&ny);

      Ex=(float **)malloc(nx*sizeof(float *));      
      for(i=0; i<nx; i++)
        Ex[i]=(float *)malloc(ny*sizeof(float));      
      dataX=(float *)malloc(nx*sizeof(float));      
      dataY=(float *)malloc(ny*sizeof(float));      

      restoreFloatArray(fileName,"X",dataX,nx);
      restoreFloatArray(fileName,"Y",dataY,ny);

      restoreField2D(fileName,dataName,Ex,nx,ny);

      sprintf(fileName,"%s_%s",argv[3],argv[4]);
      out=fopen(fileName,"w");
      for(i=0; i<nx; i+=stX)
      {
        for(j=0; j<ny; j+=stY)
          fprintf(out,"%g %g %g\n",dataX[i],dataY[j],Ex[i][j]);
        fprintf(out,"\n");
      }
      fclose(out);
      printf("%s_%s is made.\n",argv[3],argv[4]);

      for(i=0; i<nx; i++)
        free(Ex[i]);   
      free(Ex);
      free(dataX);   
      free(dataY);   
   }
}

void restoreField2D(char *fileName,char *dataName,float **data,int nx,int ny)
{
  hid_t file_id,dset_id,dataspace,memspace;
  hsize_t dimsf[2],offset[2],count[2],offset_out[2],count_out[2];
  herr_t status;
  float *field;
  int start,i,j;

  file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  dataspace=H5Dget_space(dset_id);

  //define hyperslab in dataset
  offset[0]=0;
  offset[1]=0;
  count[0]=nx;
  count[1]=ny;
  status=H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,count,NULL);

  //define memory dataspace
  dimsf[0]=nx;
  dimsf[1]=ny;
  memspace=H5Screate_simple(2,dimsf,NULL);
 
  field=(float *)malloc(nx*ny*sizeof(float));      

  //define memory hyperslab  
//  offset_out[0]=0;
//  offset_out[1]=0;
//  count_out[0]=nx;
//  count_out[1]=ny;
//  status=H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,count_out,NULL);

  //read data from hyperslab in the file into the hyperslab in memory and display 
  status=H5Dread(dset_id,H5T_NATIVE_FLOAT,memspace,dataspace,H5P_DEFAULT,field);
  start=0;
  for(i=0; i<nx; i++)
  {
    for(j=0; j<ny; j++)
      data[i][j]=field[start+j];
    start+=ny;
  }
       
  free(field);
  H5Dclose(dset_id);
  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Fclose(file_id);
}


void restoreFloatArray(char *fileName,char *dataName,float *data,int totalCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=totalCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);

  filespace=H5Screate_simple(1,metaDim,NULL);
  status=H5Dread(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}


void restoreHyperSlab2D(char *fileName,char *dataName,float *data,int nx,int ny,int pick)
{
  hid_t file_id,dset_id,dataspace,memspace;
  hsize_t dimsf[2],offset[2],count[2],offset_out[2],count_out[2];
  herr_t status;

  file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  dataspace=H5Dget_space(dset_id);

  //define hyperslab in dataset
  offset[0]=0;
  offset[1]=pick;
  count[0]=nx;
  count[1]=1;
  status=H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,count,NULL);

  //define memory dataspace
  dimsf[0]=nx;
  dimsf[1]=1;
  memspace=H5Screate_simple(2,dimsf,NULL);

  //define memory hyperslab  
//  offset_out[0]=0;
//  offset_out[1]=0;
//  count_out[0]=nx;
//  count_out[1]=ny;
//  status=H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,count_out,NULL);

  //read data from hyperslab in the file into the hyperslab in memory and display 
  status=H5Dread(dset_id,H5T_NATIVE_FLOAT,memspace,dataspace,H5P_DEFAULT,data);

  H5Dclose(dset_id);
  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Fclose(file_id);
}


void saveIntMeta(char *fileName,char *dataName,int *data)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=1;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

