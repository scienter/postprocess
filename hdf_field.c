// From x-y ascii data from xgrafix, reconstruct the density profile
#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "math.h"

void restoreFloatArray(char *fileName,char *dataName,float *data,int totalCnt);
void saveIntMeta(char *fileName,char *dataName,int *data);
void restoreHyperSlab2D(char *fileName,char *dataName,float *data,int nx,int ny,int pick);
void restoreField2D(char *fileName,char *dataName,float ***data,int nx,int ny);
void restoreField3D(char *fileName,char *dataName,float ***data,int nx,int ny,int nz,int offsetY);

main(int argc, char *argv[])
{
   FILE *in,*out;
   char fileName[100],dataName[100];
   int i,j,k,nx,ny,nz,mode,dimension,stX,stY,stZ;
   float ***field,*dataX,*dataY,*dataZ;

   if(argc < 4)
   {
      printf("hdf_field mode dimension \n");
      printf("mode:0 => fileName(.h5) dataName strideX strideY strideZ\n");
      exit(0);
   }
   mode=atoi(argv[1]);
   dimension=atoi(argv[2]);

   if(mode==0) 
   {
      sprintf(fileName,"%s.h5",argv[3]);
      sprintf(dataName,"%s",argv[4]);
      stX=atoi(argv[5]);
      stY=atoi(argv[6]);
      stZ=atoi(argv[7]);

      if(dimension>1)
      {
        saveIntMeta(fileName,"nx",&nx);
        saveIntMeta(fileName,"ny",&ny);
        nz=1;
      }
      else 	;
      if(dimension>2)
        saveIntMeta(fileName,"nz",&nz); 
      else	;

      field=(float ***)malloc(nx*sizeof(float **));      
      for(i=0; i<nx; i++)
      {
        field[i]=(float **)malloc(ny*sizeof(float *));  
        for(j=0; j<ny; j++)
          field[i][j]=(float *)malloc(nz*sizeof(float ));  
      }
      dataX=(float *)malloc(nx*sizeof(float));      
      dataY=(float *)malloc(ny*sizeof(float));      
      dataZ=(float *)malloc(nz*sizeof(float));      

      if(dimension==2)
      {        
        restoreFloatArray(fileName,"X",dataX,nx);
        restoreFloatArray(fileName,"Y",dataY,ny);

        restoreField2D(fileName,dataName,field,nx,ny);

        k=0;
        sprintf(fileName,"%s_%s",argv[3],argv[4]); 
        out=fopen(fileName,"w");
        for(i=0; i<nx; i+=stX)
        {
          for(j=0; j<ny; j+=stY)
            fprintf(out,"%g %g %g\n",dataX[i],dataY[j],field[i][j][k]);
          fprintf(out,"\n");
        }
        fclose(out);
        printf("%s_%s is made.\n",argv[3],argv[4]);
      }
      else if(dimension==3)
      {        
        restoreFloatArray(fileName,"X",dataX,nx);
        restoreFloatArray(fileName,"Y",dataY,ny);
        restoreFloatArray(fileName,"Z",dataZ,nz);

        for(j=0; j<ny; j+=stY)
        {
          restoreField3D(fileName,dataName,field,nx,ny,nz,j);
/*
          sprintf(fileName,"%s_%s_%d",argv[3],argv[4],j);
          out=fopen(fileName,"w");
          for(i=0; i<nx; i+=stX)
          {
            for(k=0; k<nz; k+=stZ)
              fprintf(out,"%g %g %g\n",dataX[i],dataZ[k],field[i][j][k]);
            fprintf(out,"\n");
          }
          fclose(out);
          printf("%s_%s is made.\n",argv[3],argv[4]);
*/
        }

      }
      else	;





      for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
          free(field[i][j]);   
      for(i=0; i<nx; i++)
        free(field[i]);   
      free(field);
      free(dataX);   
      free(dataY);   
      free(dataZ);   
   }
}

void restoreField3D(char *fileName,char *dataName,float ***data,int nx,int ny,int nz,int offsetY)
{
  hid_t file_id,dset_id,dataspace,memspace;
  hsize_t dimsf[3],offset[3],count[3],offset_out[3],count_out[3];
  herr_t status;
  float *field;
  int start,i,j,k;

  file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);

  //define memory dataspace
  dimsf[1]=nx;
  dimsf[0]=ny;
  dimsf[2]=nz;
  filespace=H5Screate_simple(3,dimsf,NULL);
 
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  dataspace=H5Dget_space(dset_id);

  //define hyperslab in dataset
  offset[1]=0;		//x
  offset[0]=offsetY;	//y
  offset[2]=0;		//z
  count[1]=nx;
  count[0]=1;
  count[2]=nz;
  memspace=H5Screate_simple(3,count,NULL);
  status=H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,count,NULL);

  field=(float *)malloc(nx*1*nz*sizeof(float));      

  //define memory hyperslab  
//  offset_out[0]=0;
//  offset_out[1]=0;
//  count_out[0]=nx;
//  count_out[1]=ny;
//  status=H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,count_out,NULL);

  //read data from hyperslab in the file into the hyperslab in memory and display 
  status=H5Dread(dset_id,H5T_NATIVE_FLOAT,memspace,dataspace,H5P_DEFAULT,field);
  start=0;
  j=offsetY;
  for(i=0; i<nx; i++)
  {
    for(k=0; k<nz; k++)
      data[i][j][k]=field[start+k];
    start+=nz;
  }
       
  free(field);
  H5Dclose(dset_id);
  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Fclose(file_id);
}

void restoreField2D(char *fileName,char *dataName,float ***data,int nx,int ny)
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
      data[i][j][0]=field[start+j];
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

