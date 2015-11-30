#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "malloc.h"
#include "hdf5.h"

main(int argc, char *argv[])
{
   float x,y,z,ux,uy,uz,id,gamma,core,minX,maxX,maxN;
   float minY,minZ,maxY,maxZ,xprime,yprime,zprime,minPx,maxPx;
   float kp,v0,gamma0,dt,bias,dx,dy,dz,rangeX,rangeY,rangeZ;
   float coefX[2],coefY[2],coefZ[2];
   float ***n;
   float *xx,*yy,*zz;
   FILE *in;
   int i,j,k,numX,numY,numZ,step,indexX,indexY,indexZ,cnt;
   int maxI,maxJ,maxK,rnk,cores;
   char name[100];


   if(argc < 12)   { 
       printf("density step numX numY numZ minX maxX minY maxY minZ maxZ minPx maxPx\n");
       exit(0);  
   }

   step=atoi(argv[1]);
   numX=atoi(argv[2]);
   numY=atoi(argv[3]);
   numZ=atoi(argv[4]);
   minX=atof(argv[5]);
   maxX=atof(argv[6]);
   minY=atof(argv[7]);
   maxY=atof(argv[8]);
   minZ=atof(argv[9]);
   maxZ=atof(argv[10]);
   minPx=atof(argv[11]);
   maxPx=atof(argv[12]);

   kp=4.21e5;
   gamma0=14.931;
//   dt=1.67e-16;
   dt=8.33e-17;
   v0=1.0-0.5/gamma0/gamma0;
   bias=step*dt*v0*kp*3e8;
//   minX=minX*kp-bias;
//   maxX=maxX*kp-bias;
//   minY=minY*kp;
//   maxY=maxY*kp;
//   minZ=minZ*kp;
//   maxZ=maxZ*kp;

   rangeX=maxX-minX;
   rangeY=maxY-minY;
   rangeZ=maxZ-minZ;
   dx=(maxX-minX)/(float)(numX);
   dy=(maxY-minY)/(float)(numY);
   dz=(maxZ-minZ)/(float)(numZ);

   n=(float ***)malloc((numX+1)*sizeof(float **));
   for(i=0; i<numX+1; i++)
   {
     n[i]=(float **)malloc((numY+1)*sizeof(float *));
     for(j=0; j<numY+1; j++)
      n[i][j]=(float *)malloc((numZ+1)*sizeof(float ));
   }
   xx=(float *)malloc((numX+1)*sizeof(float ));
   yy=(float *)malloc((numY+1)*sizeof(float ));
   zz=(float *)malloc((numZ+1)*sizeof(float ));

   for(i=0; i<numX+1; i++)
     for(j=0; j<numY+1; j++)
       for(k=0; k<numZ+1; k++)
         n[i][j][k]=0.0;
   for(i=0; i<numX+1; i++)
     xx[i]=dx*i+minX;
   for(j=0; j<numY+1; j++)
     yy[j]=dy*j+minY;
   for(k=0; k<numZ+1; k++)
     zz[k]=dz*k+minZ;

     sprintf(name,"0Particle%d",step);
     in=fopen(name,"r");
     while(fscanf(in,"%g %g %g %g %g %g %g %g %g"
                  ,&x,&y,&z,&ux,&uy,&uz,&gamma,&id,&core)!=EOF) 
     {
       xprime=x*kp-bias;
       yprime=y*kp;
       zprime=z*kp;
       xprime=x;
       yprime=y;
       zprime=z;
       indexX=(int)((xprime-minX)/dx);
       indexY=(int)((yprime-minY)/dy);
       indexZ=(int)((zprime-minZ)/dz);
       coefX[1]=(xprime-minX)/dx-indexX;
       coefX[0]=1.0-coefX[1];
       coefY[1]=(yprime-minY)/dy-indexY;
       coefY[0]=1.0-coefY[1];
       coefZ[1]=(zprime-minZ)/dz-indexZ;
       coefZ[0]=1.0-coefZ[1];
       if(indexX>=0 && indexX<numX &&
          indexY>=0 && indexY<numY &&
          indexZ>=0 && indexZ<numZ &&
          ux>minPx  && ux<maxPx)
       {
         for(i=0; i<2; i++)
           for(j=0; j<2; j++)
             for(k=0; k<2; k++)
               n[indexX+i][indexY+j][indexZ+k]+=coefX[i]*coefY[j]*coefZ[k];
       }
     }
     fclose(in);

    float *field;
    int start;
    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[3],count[3],offset[3],dimy[1],dimx[1],dimz[1];
   
    plist_id=H5Pcreate(H5P_FILE_ACCESS);
//    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
    sprintf(name,"density%d.h5",step);
    file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
    H5Pclose(plist_id);

    dimsf[0]=numY+1;      
    dimsf[1]=numX+1;     
    dimsf[2]=numZ+1;     
    filespace=H5Screate_simple(3,dimsf,NULL);

    field = (float *)malloc((numX+1)*(numY+1)*(numZ+1)*sizeof(float ));

    dset_id=H5Dcreate2(file_id,"density",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
//    subfilespace=H5Dget_space(dset_id);
//    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
    start=0;
    for(j=0; j<numY+1; j++)
      for(i=0; i<numX+1; i++)
      {
        for(k=0; k<numZ+1; k++)
          field[start+k]=n[i][j][k];
        start+=numZ+1;
      }
//    plist_id=H5Pcreate(H5P_DATASET_XFER);
//    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,field);
//    H5Pclose(plist_id);
//    H5Sclose(subfilespace);
    H5Dclose(dset_id);
    H5Sclose(filespace);

    dimx[0]=numX+1;
    dimy[0]=numY+1;
    dimz[0]=numZ+1;
    filespace=H5Screate_simple(1,dimy,NULL);
    dset_id=H5Dcreate2(file_id,"/Y",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,yy);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    filespace=H5Screate_simple(1,dimx,NULL);
    dset_id=H5Dcreate2(file_id,"/X",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,xx);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    filespace=H5Screate_simple(1,dimz,NULL);
    dset_id=H5Dcreate2(file_id,"/Z",H5T_NATIVE_FLOAT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,zz);
    H5Dclose(dset_id);
    H5Sclose(filespace);

    H5Fclose(file_id);
    printf("density%d.h is made.\n",step);

/*
   sprintf(name,"histoYZ%d",step);
   in=fopen(name,"w");
   for(j=0; j<numY+1; j++)
   {
     for(k=0; k<numZ+1; k++)
       fprintf(in,"%g %g %g\n",yy[j],zz[k],n[maxI][j][k]);
     fprintf(in,"\n");
   }
   fclose(in);
   printf("%s is made.\n",name);

   sprintf(name,"histoXY%d",step);
   in=fopen(name,"w");
   for(i=0; i<numX+1; i++)
   {
     for(j=0; j<numY+1; j++)
       fprintf(in,"%g %g %g\n",xx[i],yy[j],n[i][j][maxK]);
     fprintf(in,"\n");
   }
   fclose(in);
   printf("%s is made.\n",name);

   sprintf(name,"histoXZ%d",step);
   in=fopen(name,"w");
   for(i=0; i<numX+1; i++)
   {
     for(k=0; k<numZ+1; k++)
       fprintf(in,"%g %g %g\n",xx[i],zz[k],n[i][maxJ][k]);
     fprintf(in,"\n");
   }
   fclose(in);
   printf("%s is made.\n",name);
*/


   FILE *xmf = 0;
   sprintf(name,"density%d.xmf",step);
   xmf = fopen(name,"w");
   fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
   fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
   fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
   fprintf(xmf, " <Domain>\n");
   fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
   sprintf(name,"density");
   //3D
   fprintf(xmf, "     <Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"%d %d %d\"/>\n",numY+1,numX+1,numZ+1);
   fprintf(xmf, "     <Geometry GeometryType=\"VXVYVZ\">\n");
   fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",numZ+1);
   fprintf(xmf, "        %s%d.h5:/Z\n",name,step);
   fprintf(xmf, "       </DataItem>\n");
   fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",numX+1);
   fprintf(xmf, "        %s%d.h5:/X\n",name,step);
   fprintf(xmf, "       </DataItem>\n");
   fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",numY+1);
   fprintf(xmf, "        %s%d.h5:/Y\n",name,step);
   fprintf(xmf, "       </DataItem>\n");
   fprintf(xmf, "     </Geometry>\n");   
   fprintf(xmf, "     <Attribute Name=\"density\" AttributeType=\"Scalar\" Center=\"Node\">\n");
   fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",numZ+1,numX+1,numY+1);
   fprintf(xmf, "        %s%d.h5:/density\n",name,step);
   fprintf(xmf, "       </DataItem>\n");
   fprintf(xmf, "     </Attribute>\n");
   fprintf(xmf, "   </Grid>\n");
   fprintf(xmf, " </Domain>\n");
   fprintf(xmf, "</Xdmf>\n");
   fclose(xmf);
   printf("density%d.xmf is made.\n",step);

   free(xx);
   free(yy);
   free(zz);

   for(i=0; i<numX+1; i++)
   {
     for(j=0; j<numY+1; j++)
      free(n[i][j]);
     free(n[i]);
   }
   free(n);
   free(field);

}
