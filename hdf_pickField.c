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

main(int argc, char *argv[])
{
   FILE *in,*out;
   char fileName[100],dataName[100];
   int i,nx,ny,pick,mode,dimension;
   float pickY;
   float *Ex,*Ey,*Ez,*dataX,*dataY;

   if(argc < 4)
   {
      printf("hdf_pickField mode dimension \n");
      printf("mode:0 => fileName(.h5) pickY\n");
      exit(0);
   }
   mode=atoi(argv[1]);
   dimension=atoi(argv[2]);

   if(mode==0 && dimension==2) 
   {
      sprintf(fileName,"%s.h5",argv[3]);
      pickY=atof(argv[4]);

      saveIntMeta(fileName,"nx",&nx);
      saveIntMeta(fileName,"ny",&ny);

      Ex=(float *)malloc(nx*sizeof(float));      
      Ey=(float *)malloc(nx*sizeof(float));      
      Ez=(float *)malloc(nx*sizeof(float));      
      dataX=(float *)malloc(nx*sizeof(float));      
      dataY=(float *)malloc(ny*sizeof(float));      

      restoreFloatArray(fileName,"X",dataX,nx);
      restoreFloatArray(fileName,"Y",dataY,ny);

      pick=0;
      i=0;
      while(pick==0)
      {
        if(dataY[i]>=pickY && pickY<dataY[i+1])
          pick=i;
        else	;
        i++;
      }
      if(pick==0)    
        printf("pickY is not in the region!,data\n");
      else 	;
      restoreHyperSlab2D(fileName,"Ex",Ex,nx,ny,pick);
      restoreHyperSlab2D(fileName,"Ey",Ey,nx,ny,pick);
      restoreHyperSlab2D(fileName,"Ez",Ez,nx,ny,pick);

      sprintf(fileName,"cen%s",argv[3]);
      out=fopen(fileName,"w");
      for(i=0; i<nx; i++)
        fprintf(out,"%g %g %g %g\n",dataX[i],Ex[i],Ey[i],Ez[i]);
      fclose(out);
      printf("%s is made.\n",fileName);

      free(Ex);   
      free(Ey);   
      free(Ez);   
      free(dataX);   
      free(dataY);   
   }
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

void saveXmf(char *fileName,int numX,int numY,int numZ)
{
   FILE *xmf = 0;
   char name[100];

   sprintf(name,"%s.xmf",fileName);
   xmf = fopen(name,"w");
   fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
   fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
   fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
   fprintf(xmf, " <Domain>\n");
   fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
   sprintf(name,"density");
   //3D
   fprintf(xmf, "     <Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"%d %d %d\"/>\n",numY,numX,numZ);
   fprintf(xmf, "     <Geometry GeometryType=\"VXVYVZ\">\n");
   fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",numZ);
   fprintf(xmf, "        %s:/Z\n",fileName);
   fprintf(xmf, "       </DataItem>\n");
   fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",numX);
   fprintf(xmf, "        %s:/X\n",fileName);
   fprintf(xmf, "       </DataItem>\n");
   fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",numY);
   fprintf(xmf, "        %s:/Y\n",fileName);
   fprintf(xmf, "       </DataItem>\n");
   fprintf(xmf, "     </Geometry>\n");   
   fprintf(xmf, "     <Attribute Name=\"density\" AttributeType=\"Scalar\" Center=\"Node\">\n");
   fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",numZ,numX,numY);
   fprintf(xmf, "        %s:/density\n",fileName);
   fprintf(xmf, "       </DataItem>\n");
   fprintf(xmf, "     </Attribute>\n");
   fprintf(xmf, "   </Grid>\n");
   fprintf(xmf, " </Domain>\n");
   fprintf(xmf, "</Xdmf>\n");
   fclose(xmf);
}

void saveDataHDF3D(char *fileName,float ***ne,int numX,int numY,int numZ,float *xc,float *yc,float *zc)
{
    float *field;
    int i,j,k,start;
    float xx,yy,zz;
    char name[100];

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[3],count[3],offset[3],dimy[1],dimx[1],dimz[1];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    sprintf(name,"den%s",fileName);
    file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);

    field=(float *)malloc((numX+1)*(numY+1)*(numZ+1)*sizeof(float));
    dimsf[0]=numY+1;
    dimsf[1]=numX+1;
    dimsf[2]=numZ+1;
    filespace=H5Screate_simple(3,dimsf,NULL);
    dset_id=H5Dcreate2(file_id,"density",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    start=0;
    for(j=0; j<numY+1; j++)
      for(i=0; i<numX+1; i++)
      {
        for(k=0; k<numZ+1; k++)
          field[start+k]=ne[i][j][k];
        start+=numZ+1;
      }
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,field);
    free(field);
    H5Dclose(dset_id);
    H5Sclose(filespace);

    //coordinate
    dimsf[0]=numX+1;
    filespace=H5Screate_simple(1,dimsf,NULL);
    dset_id=H5Dcreate2(file_id,"x",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,xc);
    H5Dclose(dset_id);
    H5Sclose(filespace);

    dimsf[0]=numY+1;
    filespace=H5Screate_simple(1,dimsf,NULL);
    dset_id=H5Dcreate2(file_id,"y",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,yc);
    H5Dclose(dset_id);
    H5Sclose(filespace);

    dimsf[0]=numZ+1;
    filespace=H5Screate_simple(1,dimsf,NULL);
    dset_id=H5Dcreate2(file_id,"z",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,zc);
    H5Dclose(dset_id);
    H5Sclose(filespace);

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

