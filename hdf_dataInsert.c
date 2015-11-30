// From x-y ascii data from xgrafix, reconstruct the density profile
#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "hdf5.h"
#include "hdf5_hl.h"



void main(int argc, char *argv[])
{
   char *file,*name;
   char dataName[100],fileName[100];
   int data;

   void saveIntMeta(char *fileName,char *dataName,int *data);

   if(argc < 3)
   {
      printf("hdf_dataInsert file dataName data\n");
      exit(0);
   }
   file=argv[1];
   name=argv[2];
   data=atoi(argv[3]);

printf("%s, %s\n",file,name);
   sprintf(fileName,"%s.h5",file);
   sprintf(dataName,"%s",name);
printf("%s, %s\n",fileName,dataName);
   saveIntMeta(fileName,dataName,&data);

}

void saveIntMeta(char *fileName,char *dataName,int *data)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=1;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}


