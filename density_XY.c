#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "malloc.h"

main(int argc, char *argv[])
{
   float x,y,z,ux,uy,uz,id,gamma,core,minX,maxX,maxN,den;
   float minY,minZ,maxY,maxZ,xprime,yprime,zprime;
   float kp,v0,gamma0,dt,bias,dx,dy,dz,rangeX,rangeY,rangeZ;
   float coefX[2],coefY[2],coefZ[2];
   float **n;
   float *xx,*yy,*zz;
   FILE *in,*out;
   int i,j,k,numX,numY,numZ,step,indexX,indexY,indexZ,cnt,mode;
   int maxI,maxJ,maxK,rnk,cores;
   char name[100];


   if(argc < 12)   { 
       printf("density_XY mode step cores numX numY numZ minX maxX minY maxY minZ maxZ\n");
       printf("mode 12 : saving with x and y axis\n");
       printf("mode 13 : saving with x and z axis\n");
       exit(0);  
   }

   mode=atoi(argv[1]);
   step=atoi(argv[2]);
   cores=atoi(argv[3]);
   numX=atoi(argv[4]);
   numY=atoi(argv[5]);
   numZ=atoi(argv[6]);
   minX=atof(argv[7]);
   maxX=atof(argv[8]);
   minY=atof(argv[9]);
   maxY=atof(argv[10]);
   minZ=atof(argv[11]);
   maxZ=atof(argv[12]);

//   kp=4.21e5;
//   gamma0=14.931;
//   dt=1.67e-16;
//   v0=1.0-0.5/gamma0/gamma0;
//   bias=step*dt*v0*kp*3e8;
//   minX=minX*kp-bias;
//   maxX=maxX*kp-bias;
//   minY=minY*kp;
//   maxY=maxY*kp;
//   minZ=minZ*kp;
//   maxZ=maxZ*kp;

   switch (mode) {
   case 12 : 	 	//save with x and y axes
     rangeX=maxX-minX;
     rangeY=maxY-minY;
     dx=(maxX-minX)/(float)(numX);
     dy=(maxY-minY)/(float)(numY);

     n=(float **)malloc((numX+1)*sizeof(float *));
     for(i=0; i<numX+1; i++)
       n[i]=(float *)malloc((numY+1)*sizeof(float ));

     xx=(float *)malloc((numX+1)*sizeof(float ));
     yy=(float *)malloc((numY+1)*sizeof(float ));

     for(i=0; i<numX+1; i++)
       for(j=0; j<numY+1; j++)
         n[i][j]=0.0;
     for(i=0; i<numX+1; i++)
       xx[i]=dx*i+minX;
     for(j=0; j<numY+1; j++)
       yy[j]=dy*j+minY;

     for(rnk=0; rnk<cores; rnk++)
     {
       sprintf(name,"0Particle%d_%d",step,rnk);
       in=fopen(name,"r");
       while(fscanf(in,"%g %g %g %g %g %g %g %g %g"
                  ,&x,&y,&z,&ux,&uy,&uz,&gamma,&id,&core)!=EOF) 
       {
//         xprime=x*kp-bias;
//         yprime=y*kp;
         indexX=(int)((x-minX)/dx);
         indexY=(int)((y-minY)/dy);
         coefX[1]=(x-minX)/dx-indexX;
         coefX[0]=1.0-coefX[1];
         coefY[1]=(y-minY)/dy-indexY;
         coefY[0]=1.0-coefY[1];
         if(indexX>=0 && indexX<numX &&
            indexY>=0 && indexY<numY &&
            z>=minZ && z<maxZ)
         {
           for(i=0; i<2; i++)
             for(j=0; j<2; j++)
               n[indexX+i][indexY+j]+=coefX[i]*coefY[j];
         }
       }
       fclose(in);
       printf("core%d is done.\n",rnk);
     }

     sprintf(name,"denXY%d",step);
     out=fopen(name,"w");
     for(i=0; i<numX+1; i++)
     {
       for(j=0; j<numY+1; j++)
       {
         x=xx[i];
         y=yy[j];
         den=n[i][j];
         fprintf(out,"%g %g %g\n",x,y,den);
       }
       fprintf(out,"\n");
     }
     fclose(out);
     printf("denXY%d is made.\n",step);

     free(xx);
     free(yy);
     for(i=0; i<numX+1; i++)
       free(n[i]);
     free(n);
     break;
   case 13 : 	 	//save with x and z axes
     rangeX=maxX-minX;
     rangeZ=maxZ-minZ;
     dx=(maxX-minX)/(float)(numX);
     dz=(maxZ-minZ)/(float)(numZ);

     n=(float **)malloc((numX+1)*sizeof(float *));
     for(i=0; i<numX+1; i++)
       n[i]=(float *)malloc((numZ+1)*sizeof(float ));

     xx=(float *)malloc((numX+1)*sizeof(float ));
     zz=(float *)malloc((numZ+1)*sizeof(float ));

     for(i=0; i<numX+1; i++)
       for(k=0; k<numZ+1; k++)
         n[i][k]=0.0;
     for(i=0; i<numX+1; i++)
       xx[i]=dx*i+minX;
     for(k=0; k<numZ+1; k++)
       zz[k]=dz*k+minZ;

     for(rnk=0; rnk<cores; rnk++)
     {
       sprintf(name,"0Particle%d_%d",step,rnk);
       in=fopen(name,"r");
       while(fscanf(in,"%g %g %g %g %g %g %g %g %g"
                  ,&x,&y,&z,&ux,&uy,&uz,&gamma,&id,&core)!=EOF) 
       {
//         xprime=x*kp-bias;
//         zprime=z*kp;
         indexX=(int)((x-minX)/dx);
         indexZ=(int)((z-minZ)/dz);
         coefX[1]=(x-minX)/dx-indexX;
         coefX[0]=1.0-coefX[1];
         coefZ[1]=(z-minZ)/dz-indexZ;
         coefZ[0]=1.0-coefZ[1];
         if(indexX>=0 && indexX<numX &&
            indexZ>=0 && indexZ<numZ &&
            y>=minY && y<maxY)
         {
           for(i=0; i<2; i++)
             for(k=0; k<2; k++)
               n[indexX+i][indexZ+k]+=coefX[i]*coefZ[k];
         }
       }
       fclose(in);
       printf("core%d is done.\n",rnk);
     }

     sprintf(name,"denXZ%d",step);
     out=fopen(name,"w");
     for(i=0; i<numX+1; i++)
     {
       for(k=0; k<numZ+1; k++)
       {
         x=xx[i];
         z=zz[k];
         den=n[i][k];
         fprintf(out,"%g %g %g\n",x,z,den);
       }
       fprintf(out,"\n");
     }
     fclose(out);
     printf("denXZ%d is made.\n",step);

     free(xx);
     free(zz);
     for(i=0; i<numX+1; i++)
       free(n[i]);
     free(n);
     break;
   }		//End of switch (mode)
 

}
