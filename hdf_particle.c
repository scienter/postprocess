// From x-y ascii data from xgrafix, reconstruct the density profile
#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "hdf5.h"
#include "hdf5_hl.h"


typedef struct _Param  {
  struct _Ptcl *head;
} Param;


typedef struct _Ptcl  {
  int id;
  int core; 
  struct _Ptcl *next;
} Ptcl;


void main(int argc, char *argv[])
{
   FILE *out,*in;
   char dataName[100],fileName[100],name[100];
   int i,ii,step,totalCnt,targetStep,dimension,mode,species,cnt,testCnt;
   int index,initial,final,timeStep,id,core,minCore,maxCore,flag=0;
   float x,y,z,px,py,pz,minX,maxX,rangeX,minPx,incline,minY,maxY;
   float *dataX,*dataY,*dataZ,*dataPx,*dataPy,*dataPz;
   Param D;
   Ptcl *p;

   int *selectIndex,*selectCore,*dataIndex,*dataCore;
   void restoreIntMeta(char *fileName,char *dataName,int *data);
   void restoreIntArray(char *fileName,char *dataName,int *data,int totalCnt);
   void restoreFloatArray(char *fileName,char *dataName,float *data,int totalCnt);
   void createParticle(Param *D,int *selectIndex,int *selectCore,int cnt);
   int testingData(Param *D,int id,int core);
   void removeData(Param *D);

   if(argc < 6)
   {
      printf("hdf_particle mode dimension species\n");
      printf("mode(1) for hdf file : targetStep minX maxX minY maxY\n");
      printf("mode(2) for hdf file : targetStep incline rangeX minPx\n");
      printf("mode(3) for txt file : fileName minX maxX minY maxY\n");
      printf("mode(4) for txt file : fileName incline rangeX minPx\n");
      printf("mode(5) using 'idsample' file : initial final timeStep\n");
      exit(0);
   }
   mode=atoi(argv[1]);
   dimension=atoi(argv[2]);
   species=atoi(argv[3]);

   if(dimension==2 && mode==1)  //2D
   {
     targetStep=atoi(argv[4]);
     minX=atof(argv[5]);
     maxX=atof(argv[6]);
     minY=atof(argv[7]);
     maxY=atof(argv[8]);

     sprintf(fileName,"%dParticle%d.h5",species,targetStep);
     sprintf(dataName,"/totalCnt");
     restoreIntMeta(fileName,dataName,&totalCnt);

     dataX=(float *)malloc(totalCnt*sizeof(float));
     dataY=(float *)malloc(totalCnt*sizeof(float));
     dataPx=(float *)malloc(totalCnt*sizeof(float));
     dataPy=(float *)malloc(totalCnt*sizeof(float));
     dataPz=(float *)malloc(totalCnt*sizeof(float));
     dataIndex=(int *)malloc(totalCnt*sizeof(int));
     dataCore=(int *)malloc(totalCnt*sizeof(int));

     sprintf(dataName,"/x");
     restoreFloatArray(fileName,dataName,dataX,totalCnt);
     sprintf(dataName,"/y");
     restoreFloatArray(fileName,dataName,dataY,totalCnt);
     sprintf(dataName,"/px");
     restoreFloatArray(fileName,dataName,dataPx,totalCnt);
     sprintf(dataName,"/py");
     restoreFloatArray(fileName,dataName,dataPy,totalCnt);
     sprintf(dataName,"/pz");
     restoreFloatArray(fileName,dataName,dataPz,totalCnt);
     sprintf(dataName,"/index");
     restoreIntArray(fileName,dataName,dataIndex,totalCnt);
     sprintf(dataName,"/core");
     restoreIntArray(fileName,dataName,dataCore,totalCnt);

     sprintf(name,"%dParticle%d",species,targetStep);
     out=fopen(name,"w");
     cnt=0;
     for (i=0; i<totalCnt; i++)
     {
       x=dataX[i];
       y=dataY[i];
       px=dataPx[i];
       py=dataPy[i];
       pz=dataPz[i];
       id=dataIndex[i];
       core=dataCore[i];
       if(x>minX && x< maxX && y>minY && y<maxY)
       {
         cnt++;
         fprintf(out,"%g %g %g %g %g %g %d %d\n",x,y,0.0,px,py,pz,id,core);
       }
     }
     fclose(out);
     printf("%dParticle%d is made.\n",species,targetStep);

     free(dataX);
     free(dataY);
     free(dataPx);
     free(dataPy);
     free(dataPz);
     free(dataIndex);
     free(dataCore);
   }
 
   else if(dimension==2 && mode==2)  //2D
   {
     targetStep=atoi(argv[4]);
     incline=atof(argv[5]);
     maxX=incline*targetStep;
     rangeX=atof(argv[6]);
     minX=maxX-rangeX;
     minPx=atof(argv[7]);

     sprintf(fileName,"%dParticle%d.h5",species,targetStep);
     sprintf(dataName,"/totalCnt");
     restoreIntMeta(fileName,dataName,&totalCnt);

     dataX=(float *)malloc(totalCnt*sizeof(float));
     dataY=(float *)malloc(totalCnt*sizeof(float));
     dataPx=(float *)malloc(totalCnt*sizeof(float));
     dataPy=(float *)malloc(totalCnt*sizeof(float));
     dataPz=(float *)malloc(totalCnt*sizeof(float));
     dataIndex=(int *)malloc(totalCnt*sizeof(int));
     dataCore=(int *)malloc(totalCnt*sizeof(int));

     sprintf(dataName,"/x");
     restoreFloatArray(fileName,dataName,dataX,totalCnt);
     sprintf(dataName,"/y");
     restoreFloatArray(fileName,dataName,dataY,totalCnt);
     sprintf(dataName,"/px");
     restoreFloatArray(fileName,dataName,dataPx,totalCnt);
     sprintf(dataName,"/py");
     restoreFloatArray(fileName,dataName,dataPy,totalCnt);
     sprintf(dataName,"/pz");
     restoreFloatArray(fileName,dataName,dataPz,totalCnt);
     sprintf(dataName,"/index");
     restoreIntArray(fileName,dataName,dataIndex,totalCnt);
     sprintf(dataName,"/core");
     restoreIntArray(fileName,dataName,dataCore,totalCnt);

     sprintf(name,"%dParticle%d",species,targetStep);
     out=fopen(name,"w");
     cnt=0;
     for (i=0; i<totalCnt; i++)
     {
       x=dataX[i];
       y=dataY[i];
       px=dataPx[i];
       py=dataPy[i];
       pz=dataPz[i];
       id=dataIndex[i];
       core=dataCore[i];
       if(x>minX && px>minPx)
       {
         cnt++;
         fprintf(out,"%g %g %g %g %g %g %d %d\n",x,y,0.0,px,py,pz,id,core);
       }
     }
     fclose(out);
     printf("%dParticle%d is made. \n",species,targetStep);

     free(dataX);
     free(dataY);
     free(dataPx);
     free(dataPy);
     free(dataPz);
     free(dataIndex);
     free(dataCore);
   }
   else if(dimension==3 && mode==2)  //3D
   {
     targetStep=atoi(argv[4]);
     incline=atof(argv[5]);
     maxX=incline*targetStep;
     rangeX=atof(argv[6]);
     minX=maxX-rangeX;
     minPx=atof(argv[7]);

     sprintf(fileName,"%dParticle%d.h5",species,targetStep);
     sprintf(dataName,"/totalCnt");
     restoreIntMeta(fileName,dataName,&totalCnt);

     dataX=(float *)malloc(totalCnt*sizeof(float));
     dataY=(float *)malloc(totalCnt*sizeof(float));
     dataZ=(float *)malloc(totalCnt*sizeof(float));
     dataPx=(float *)malloc(totalCnt*sizeof(float));
     dataPy=(float *)malloc(totalCnt*sizeof(float));
     dataPz=(float *)malloc(totalCnt*sizeof(float));
     dataIndex=(int *)malloc(totalCnt*sizeof(int));
     dataCore=(int *)malloc(totalCnt*sizeof(int));

     sprintf(dataName,"/x");
     restoreFloatArray(fileName,dataName,dataX,totalCnt);
     sprintf(dataName,"/y");
     restoreFloatArray(fileName,dataName,dataY,totalCnt);
     sprintf(dataName,"/z");
     restoreFloatArray(fileName,dataName,dataZ,totalCnt);
     sprintf(dataName,"/px");
     restoreFloatArray(fileName,dataName,dataPx,totalCnt);
     sprintf(dataName,"/py");
     restoreFloatArray(fileName,dataName,dataPy,totalCnt);
     sprintf(dataName,"/pz");
     restoreFloatArray(fileName,dataName,dataPz,totalCnt);
     sprintf(dataName,"/index");
     restoreIntArray(fileName,dataName,dataIndex,totalCnt);
     sprintf(dataName,"/core");
     restoreIntArray(fileName,dataName,dataCore,totalCnt);

     sprintf(name,"%dParticle%d",species,targetStep);
     out=fopen(name,"w");
     cnt=0;
     for (i=0; i<totalCnt; i++)
     {
       x=dataX[i];
       y=dataY[i];
       z=dataZ[i];
       px=dataPx[i];
       py=dataPy[i];
       pz=dataPz[i];
       id=dataIndex[i];
       core=dataCore[i];
       if(x>minX && px>minPx)
       {
         cnt++;
         fprintf(out,"%g %g %g %g %g %g %d %d\n",x,y,z,px,py,pz,id,core);
       }
     }
     fclose(out);
     printf("%dParticle%d is made. \n",species,targetStep);

     free(dataX);
     free(dataY);
     free(dataZ);
     free(dataPx);
     free(dataPy);
     free(dataPz);
     free(dataIndex);
     free(dataCore);
   }

   else if(dimension==2 && mode==3)  //2D
   {
     minX=atof(argv[5]);
     maxX=atof(argv[6]);
     minY=atof(argv[7]);
     maxY=atof(argv[8]);

     cnt=0;
     in=fopen(argv[4],"r");
     out=fopen("idsample","w");
     while(fscanf(in,"%g %g %g %g %g %g %d %d",&x,&y,&z,&px,&py,&pz,&id,&core)!=EOF)
     {
       if(x>minX && x< maxX && y>minY && y<maxY)
       {
         fprintf(out,"%g %g %g %g %g %g %d %d\n",x,y,z,px,py,pz,id,core);
       }
     }
     fclose(out);
     fclose(in);
     printf("idsample is made.\n");
   }

   else if(dimension==2 && mode==4)  //2D
   {
     incline=atof(argv[5]);
     maxX=incline*targetStep;
     rangeX=atof(argv[6]);
     minX=maxX-rangeX;
     minPx=atof(argv[7]);

     cnt=0;
     in=fopen(argv[4],"r");
     out=fopen("idsample","w");
     while(fscanf(in,"%g %g %g %g %d %d",&x,&y,&z,&px,&id,&core)!=EOF)
     {
       if(x>minX && px>minPx)
       {
         fprintf(out,"%g %g %g %g %d %d\n",x,y,z,px,id,core);
       }
     }
     fclose(out);
     fclose(in);
   }

   else if(dimension==2 && mode==5)  //2D
   {
     initial=atoi(argv[4]);
     final=atoi(argv[5]);
     timeStep=atoi(argv[6]);

     cnt=0;
     in=fopen("idsample","r");
     while(fscanf(in,"%g %g %g %g %g %g %d %d",&x,&y,&z,&px,&py,&pz,&id,&core)!=EOF)
     {
       cnt++;
     }
     fclose(in);

     selectIndex=(int *)malloc(cnt*sizeof(int));
     selectCore=(int *)malloc(cnt*sizeof(int));

     minCore=200000;
     maxCore=-1;
     in=fopen("idsample","r");
     for(i=0; i<cnt; i++)
     {
       fscanf(in,"%g %g %g %g %g %g %d %d",&x,&y,&z,&px,&py,&pz,&selectIndex[i],&selectCore[i]);
       core=selectCore[i];
       if(minCore>core) minCore=core;
       else if(maxCore<core) maxCore=core;
     }
     fclose(in);
     printf("minCore=%d, maxCore=%d\n",minCore,maxCore);

     D.head=NULL; 
     for(step=initial; step<=final; step+=timeStep)
     {
       sprintf(fileName,"%dParticle%d.h5",species,step);
       sprintf(dataName,"/totalCnt");
       restoreIntMeta(fileName,dataName,&totalCnt);
       
       createParticle(&D,selectIndex,selectCore,cnt);     

       dataX=(float *)malloc(totalCnt*sizeof(float));
       dataY=(float *)malloc(totalCnt*sizeof(float));
       dataPx=(float *)malloc(totalCnt*sizeof(float));
       dataPy=(float *)malloc(totalCnt*sizeof(float));
       dataPz=(float *)malloc(totalCnt*sizeof(float));
       dataIndex=(int *)malloc(totalCnt*sizeof(int));
       dataCore=(int *)malloc(totalCnt*sizeof(int));

       sprintf(dataName,"/x");
       restoreFloatArray(fileName,dataName,dataX,totalCnt);
       sprintf(dataName,"/y");
       restoreFloatArray(fileName,dataName,dataY,totalCnt);
       sprintf(dataName,"/px");
       restoreFloatArray(fileName,dataName,dataPx,totalCnt);
       sprintf(dataName,"/py");
       restoreFloatArray(fileName,dataName,dataPy,totalCnt);
       sprintf(dataName,"/pz");
       restoreFloatArray(fileName,dataName,dataPz,totalCnt);
       sprintf(dataName,"/index");
       restoreIntArray(fileName,dataName,dataIndex,totalCnt);
       sprintf(dataName,"/core");
       restoreIntArray(fileName,dataName,dataCore,totalCnt);

       sprintf(name,"%did%d",species,step);
       out=fopen(name,"w");
       testCnt=cnt;
       for(i=0; i<totalCnt; i++)
       {
         id=dataIndex[i];
         core=dataCore[i];
         x=dataX[i];
         y=dataY[i];
         px=dataPx[i];
         py=dataPy[i];
         pz=dataPz[i];
         if(core>=minCore && core<=maxCore)
         {
           flag=testingData(&D,id,core);
           if(flag==1)
             fprintf(out,"%g %g %g %g %g %g %d %d\n",x,y,0.0,px,py,pz,id,core);
         }
       }
       fclose(out);
       printf("%did%d is saved.\n",species,step);
 
       free(dataX);
       free(dataY);
       free(dataPx);
       free(dataPy);
       free(dataPz);
       free(dataIndex);
       free(dataCore);
     }

     free(selectIndex);
     free(selectCore);
   }

}

void removeData(Param *D)
{
  int cnt,flag;
  Ptcl *p;
//here
  cnt=0;
  p=D->head;
  while(p)
  {
    cnt++;
    p=p->next;
  }
  if(cnt>0) printf("number of remained data is %d.\n",cnt);

  cnt=1;
  p=D->head;
  while(p)
  {
    D->head=p->next;
    p->next=NULL;
    free(p);
    p=D->head; 
  }

}

int testingData(Param *D,int id,int core)
{
  int cnt,flag;
  Ptcl *p,*prev;
//here
  cnt=1;
  flag=0;
  p=D->head;
  while(p)
  {
    if(cnt==1)
      prev=p;

    if(id==p->id && core==p->core)
    {
      if(cnt==1)
      {
        D->head=p->next;
        p->next=NULL;
        free(p);
        p=D->head; 
        cnt=1;
      }
      else
      {
        prev->next=p->next;
        p->next=NULL;
        free(p);       
        p=prev->next;
      }
      flag=1;
    }
    else
    {
      prev=p;
      p=p->next;
      cnt++;
    }

  }

  return flag;
}
    
void createParticle(Param *D,int *selectIndex,int *selectCore,int cnt)
{
  int i;
  Ptcl *new;
  
  i=0;
  while(i<cnt)
  {
    new=(Ptcl *)malloc(sizeof(Ptcl));
    new->next=D->head;
    D->head=new;
    new->id=selectIndex[i];
    new->core=selectCore[i];
    i++;
  }
}       

void restoreFloatArray(char *fileName,char *dataName,float *data,int totalCnt)
{
  hid_t file_id,dset_id,dspace_id,mspace_id;
  hsize_t dims[1],dimsm[1],count[1],offset[1],stride[1],block[1];
  herr_t status;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);

//  offset[0]=offSet;
//  count[0]=subCnt;
//  stride[0]=1;
//  block[0]=1;
  
  dims[0]=totalCnt;
  mspace_id = H5Screate_simple(1,dims,NULL); 
//  dspace_id = H5Dget_space(dset_id);
//  status = H5Sselect_hyperslab(dset_id,H5S_SELECT_SET,offset,stride,count,block); 
  status = H5Dread(dset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

  H5Sclose(mspace_id);
  H5Dclose(dset_id);
  H5Fclose(file_id);
}

void restoreIntArray(char *fileName,char *dataName,int *data,int totalCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=totalCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
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


