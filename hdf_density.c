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

main(int argc, char *argv[])
{
   FILE *in,*out;
   char saveFile[100],fileName[100];
   float x,y,z;
   int i,j,k,nn,Nx,Ny,Nz,Npx,mode,saveMode,dimension,core,index,NE,totalCnt;
   float xmin,xmax,ymin,ymax,zmin,zmax,minPx,maxPx,dx,dy,dz,dpx,np2c;
   float minE,maxE,dE,px,py,pz,energy,minX,weightnum,xx,yy,zz;
   float ***ne,*xc,*yc,*zc,*pxc,*nE,*xE,*dataX,*dataY,*dataZ;

   if(argc < 8)
   {
      printf("hdf_density mode dimension saveMode(0:txt,1:hdf)\n");
      printf("mode:0 => fileName(txt) minX maxX Nx minY maxY Ny np2c\n");
      printf("mode:1 => fileName(txt) minX maxX Nx minPx maxPx Npx minY maxY np2c\n");
      printf("mode:2 => fileName(txt) minE maxE NE minX weightNum\n");
      printf("mode:3 => fileName(hdf) minX maxX Nx minY maxY Ny np2c\n");
      exit(0);
   }
   mode=atoi(argv[1]);
   dimension=atoi(argv[2]);
   saveMode=atoi(argv[3]);

   if(mode==0 && dimension==2) {
      xmin=atof(argv[5]);
      xmax=atof(argv[6]);
      Nx=atof(argv[7]);
      ymin=atof(argv[8]);
      ymax=atof(argv[9]);
      Ny=atof(argv[10]);
      np2c=atof(argv[11]);

      dx = (xmax-xmin)/((float)Nx);
      dy = (ymax-ymin)/((float)Ny);

      ne=(float ***)malloc((Nx+1)*sizeof(float **));
      xc=(float *)malloc((Nx+1)*sizeof(float));
      yc=(float *)malloc((Ny+1)*sizeof(float));

      for(i=0;i<=Nx;i++) 
        ne[i]=(float **)malloc((Ny+1)*sizeof(float *));
      for(i=0;i<=Nx;i++) 
        for(j=0;j<=Ny;j++) 
          ne[i][j]=(float *)malloc((Nz+1)*sizeof(float ));
      for(i=0;i<=Nx;i++)  xc[i]=xmin+i*dx;
      for(j=0;j<=Ny;j++)  yc[j]=ymin+j*dy;

      for(i=0;i<=Nx;i++)
        for(j=0;j<=Ny;j++)
          for(k=0;k<1;k++)
            ne[i][j][k]=0.0;

      in=fopen(argv[4],"r");
      while(fscanf(in,"%g %g %g %g %g %g %d %d",&x,&y,&z,&px,&py,&pz,&index,&core)!=EOF)
      {
         i=(int)((x-xmin)/dx);
         j=(int)((y-ymin)/dy);
         k=0;
         if(0 < i && i< Nx && 0 < j && j < Ny)
         {
            ne[i][j][k]+=(xc[i+1]-x)*(yc[j+1]-y)/dx/dy*np2c;
            ne[i+1][j][k]+=(x-xc[i])*(yc[j+1]-y)/dx/dy*np2c;
            ne[i][j+1][k]+=(xc[i+1]-x)*(y-yc[j])/dx/dy*np2c;
            ne[i+1][j+1][k]+=(x-xc[i])*(y-yc[j])/dx/dy*np2c;
         }
      }
      fclose(in);

      if(saveMode==0)
      {
        sprintf(saveFile,"den%s",argv[4]);
        out=fopen(saveFile,"w");
        for(i=0; i<Nx; i++)
        {
          for(j=0; j<Ny; j++)
            fprintf(out,"%g %g %g\n",xc[i],yc[j],ne[i][j][0]/dx/dy);
          fprintf(out,"\n");
        }
        fclose(out);   
        printf("den%s is made.\n",argv[4]); 
      }

      free(xc);
      free(yc);
      for(i=0; i<=Nx; i++)
        free(ne[i]);
      free(ne);

   }

   else if(mode==1 && dimension==2) {
      xmin=atof(argv[5]);
      xmax=atof(argv[6]);
      Nx=atof(argv[7]);
      minPx=atof(argv[8]);
      maxPx=atof(argv[9]);
      Npx=atof(argv[10]);
      ymin=atof(argv[11]);
      ymax=atof(argv[12]);
      np2c=atof(argv[13]);

      dx = (xmax-xmin)/((float)Nx);
      dpx = (maxPx-minPx)/((float)Npx);

      ne=(float ***)malloc((Nx+1)*sizeof(float **));
      xc=(float *)malloc((Nx+1)*sizeof(float));
      pxc=(float *)malloc((Npx+1)*sizeof(float));

      for(i=0;i<=Nx;i++) 
        ne[i]=(float **)malloc((Ny+1)*sizeof(float *));
      for(i=0;i<=Nx;i++) 
        for(j=0;j<=Ny;j++) 
          ne[i][j]=(float *)malloc((Nz+1)*sizeof(float ));
      for(i=0;i<=Nx;i++)  xc[i]=xmin+i*dx;
      for(j=0;j<=Npx;j++)  pxc[j]=minPx+j*dpx;
      for(i=0;i<=Nx;i++)
        for(j=0;j<=Ny;j++)
          for(k=0;k<1;k++)
            ne[i][j][k]=0.0;

      in=fopen(argv[4],"r");
      while(fscanf(in,"%g %g %g %g %g %g %d %d",&x,&y,&z,&px,&py,&pz,&index,&core)!=EOF)
      {
         i=(int)((x-xmin)/dx);
         j=(int)((px-minPx)/dpx);
         if(0 < i && i< Nx && 0 < j && j < Npx && y>ymin && y<ymax)
         {
            ne[i][j][0]+=(xc[i+1]-x)*(pxc[j+1]-px)/dx/dpx*np2c;
            ne[i+1][j][0]+=(x-xc[i])*(pxc[j+1]-px)/dx/dpx*np2c;
            ne[i][j+1][0]+=(xc[i+1]-x)*(px-pxc[j])/dx/dpx*np2c;
            ne[i+1][j+1][0]+=(x-xc[i])*(px-pxc[j])/dx/dpx*np2c;
         }
      }
      fclose(in);

      if(saveMode==0)
      {
        sprintf(saveFile,"denE%s",argv[4]);
        out=fopen(saveFile,"w");
        for(i=0; i<Nx; i++)
        {
          for(j=0; j<Npx; j++)
            fprintf(out,"%g %g %g\n",xc[i],pxc[j],ne[i][j][0]);
          fprintf(out,"\n");
        }
        fclose(out);  
        printf("denE%s is made.\n",argv[4]); 
      }

      free(xc);
      free(pxc);
      for(i=0; i<=Nx; i++)
        free(ne[i]);
      free(ne);

   }

   else if(mode==2 && dimension==2) {
      minE=atof(argv[5]);
      maxE=atof(argv[6]);
      NE=atof(argv[7]);
      minX=atof(argv[8]);
      weightnum=atof(argv[9]);

      dE = (maxE-minE)/((float)NE);

      nE=(float *)malloc((NE+1)*sizeof(float));
      xE=(float *)malloc((NE+1)*sizeof(float));

      for(i=0;i<=NE;i++)
      {
        nE[i]=0.0;
        xE[i]=minE+i*dE;
      }

      in=fopen(argv[4],"r");
      while(fscanf(in,"%g %g %g %g %g %g %d %d",&x,&y,&z,&px,&py,&pz,&index,&core)!=EOF)
      {
         energy=(sqrt(1.0+px*px+py*py+pz*pz)-1.0)*938.272*weightnum; //MeV
//printf("dE=%g, energy=%g\n",dE,energy);
         i=(int)((energy-minE)/dE);
         if(0 < i && i< NE && x>minX)
         {
            nE[i]+=(xE[i+1]-energy)/dE;
            nE[i+1]+=(energy-xE[i])/dE;
         }
      }
      fclose(in);

      if(saveMode==0)
      {
        sprintf(saveFile,"spec%s",argv[4]);
        out=fopen(saveFile,"w");
        for(i=0; i<NE; i++)
          fprintf(out,"%g %g\n",xE[i],nE[i]);
        fclose(out);  
        printf("spec%s is made.\n",argv[4]); 
      }

      free(xE);
      free(nE);

   }

   else if(mode==3 && dimension==3) {
      xmin=atof(argv[5]);
      xmax=atof(argv[6]);
      Nx=atof(argv[7]);
      ymin=atof(argv[8]);
      ymax=atof(argv[9]);
      Ny=atof(argv[10]);
      np2c=atof(argv[11]);

      dx = (xmax-xmin)/((float)Nx);
      dy = (ymax-ymin)/((float)Ny);
      dz = dy;
      Nz = Ny;

      ne=(float ***)malloc((Nx+1)*sizeof(float **));
      xc=(float *)malloc((Nx+1)*sizeof(float));
      yc=(float *)malloc((Ny+1)*sizeof(float));
      zc=(float *)malloc((Nz+1)*sizeof(float));
      for(i=0;i<=Nx;i++) 
        ne[i]=(float **)malloc((Ny+1)*sizeof(float *));
      for(i=0;i<=Nx;i++) 
        for(j=0;j<=Ny;j++) 
          ne[i][j]=(float *)malloc((Nz+1)*sizeof(float ));

      for(i=0;i<=Nx;i++)
        for(j=0;j<=Ny;j++)
          for(k=0;k<=Nz;k++)
            ne[i][j][k]=0.0;
      for(i=0;i<=Nx;i++)  xc[i]=xmin+i*dx;
      for(j=0;j<=Ny;j++)  yc[j]=ymin+j*dy;
      for(k=0;k<=Nz;k++)  zc[k]=zmin+k*dz;

      sprintf(saveFile,"%s",argv[4]);
      saveIntMeta(saveFile,"totalCnt",&totalCnt);
      dataX=(float *)malloc(totalCnt*sizeof(float));      
      dataY=(float *)malloc(totalCnt*sizeof(float));      
      dataZ=(float *)malloc(totalCnt*sizeof(float));      
      restoreFloatArray(saveFile,"x",dataX,totalCnt);
      restoreFloatArray(saveFile,"y",dataY,totalCnt);
      restoreFloatArray(saveFile,"z",dataZ,totalCnt);

      for(nn=0; nn<totalCnt; nn++)
      {
        xx=dataX[nn];
        yy=dataY[nn];
        zz=dataZ[nn];
        i=(int)((xx-xmin)/dx);
        j=(int)((yy-ymin)/dy);
        k=(int)((zz-zmin)/dz);
        if(i>=0 && i<Nx && j>=0 && j<Ny && k>=0 && k<Nz)
        {
          ne[i][j][k]+=(xc[i+1]-xx)*(yc[j+1]-yy)*(zc[k+1]-zz)/dx/dy/dz*np2c;
          ne[i][j][k+1]+=(xc[i+1]-xx)*(yc[j+1]-yy)*(zz-zc[k])/dx/dy/dz*np2c;
          ne[i][j+1][k]+=(xc[i+1]-xx)*(yy-yc[j])*(zc[k+1]-zz)/dx/dy/dz*np2c;
          ne[i][j+1][k+1]+=(xc[i+1]-xx)*(yy-yc[j])*(zz-zc[k])/dx/dy/dz*np2c;
          ne[i+1][j][k]+=(xx-xc[i])*(yc[j+1]-yy)*(zc[k+1]-zz)/dx/dy/dz*np2c;
          ne[i+1][j][k+1]+=(xx-xc[i])*(yc[j+1]-yy)*(zz-zc[k])/dx/dy/dz*np2c;
          ne[i+1][j+1][k]+=(xx-xc[i])*(yy-yc[j])*(zc[k+1]-zz)/dx/dy/dz*np2c;
          ne[i+1][j+1][k+1]+=(xx-xc[i])*(yy-yc[j])*(zz-zc[k])/dx/dy/dz*np2c;
        }
        else	;
      }

      saveDataHDF3D(saveFile,ne,Nx,Ny,Nz,xc,yc,zc);
      sprintf(fileName,"den%s",saveFile);
      saveXmf(fileName,Nx+1,Ny+1,Nz+1);
//here
      
      free(dataX);
      free(dataY);
      free(dataZ);
      free(xc);
      free(yc);
      free(zc);
      for(i=0;i<=Nx;i++) 
        for(j=0;j<=Ny;j++) 
          free(ne[i][j]);
      for(i=0;i<=Nx;i++) 
        free(ne[i]);
      free(ne);

   }

   

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

