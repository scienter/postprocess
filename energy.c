#include "stdio.h"
#include "stdlib.h"
#include "math.h"
main(int argc, char *argv[])
{
   float x,y,z,ux,uy,uz,id,core,gamma,minE,maxE,energy,coefE,coefX,coef;
   float minX,maxX,dt,rangeX,superP,dx;
   int numE,numX,i,j,index,indexX,indexE,step,initial,final,timeStep,mode;
   FILE *in,*out;
   float *nE,**nnE;
   char str[100];

   if(argc < 6)
   {  printf("mode\n");
      printf("mode 0 : initial final timeStep maxE(MeV) rangeX dt(1.67e-16) superP(3.13e4) (for 2D)\n");
      printf("mode 1 : initial final timeStep maxE(MeV) rangeX dt(1.67e-16) superP(3.13e4) numX (for 3D)\n");
      printf("mode 2 : file minX maxX numX maxE(MeV) step superP(3.13e4) (for 3D and given file)\n");
      exit(0);
   }
   mode=atoi(argv[1]);

   switch (mode)  {
   case 0 :
     initial=atoi(argv[2]);
     final=atoi(argv[3]);
     timeStep=atoi(argv[4]);
     maxE=atof(argv[5]);
     rangeX=atof(argv[6]);
     dt=atof(argv[7]);
     superP=atof(argv[8]);

     numE=(int)(maxE);
     nE=(float *)malloc((numE+1)*sizeof(float ));

     for(step=initial; step<=final; step+=timeStep)
     {
       maxX=step*dt*3e8;
       minX=maxX-rangeX;
       for(i=0; i<=numE; i++)
         nE[i]=0.0;

       sprintf(str,"id%d",step);
       in = fopen(str,"r");
//   fgets(str,100,in);
       while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&x,&y,&z,&ux,&uy,&uz,&gamma,&id,&core)!=EOF)
       {      
         energy=gamma*0.511;
         index=(int)(energy);
         coef=energy-index;
         if(minX<x && x<maxX)
         {
           nE[index]+=(1.0-coef);
           nE[index+1]+=coef;
         }
       }
       fclose(in);

       sprintf(str,"energy%d",step);
       out=fopen(str,"w");
       for(i=0; i<numE; i++)
       {
         fprintf(out,"%d %g\n",i,nE[i]*superP*1.6e-19);
       }
       fclose(out); 
       printf("energy%d is made.\n",step);
     }		//End of step
     free(nE);
     break;
   
   case 1 :
     initial=atoi(argv[2]);
     final=atoi(argv[3]);
     timeStep=atoi(argv[4]);
     maxE=atof(argv[5]);
     rangeX=atof(argv[6]);
     dt=atof(argv[7]);
     superP=atof(argv[8]);
     numX=atoi(argv[9]);

     numE=(int)(maxE);     

     nnE=(float **)malloc((numE+1)*sizeof(float *));
     for(i=0; i<=numE; i++)
       nnE[i]=(float *)malloc((numX+1)*sizeof(float ));


     for(step=initial; step<=final; step+=timeStep)
     {
       maxX=step*dt*3e8;
       minX=maxX-rangeX;
       dx=rangeX/((float)numX);

       for(i=0; i<numE; i++)
         for(j=0; j<numX; j++)
           nnE[i][j]=0.0;

       sprintf(str,"id%d",step);
       in = fopen(str,"r");
//   fgets(str,100,in);
       while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&x,&y,&z,&ux,&uy,&uz,&gamma,&id,&core)!=EOF)
       {      
         energy=gamma*0.511;
         indexE=(int)(energy);
         coefE=energy-indexE;
         indexX=(int)((x-minX)/dx);
         coefX=(x-minX)/dx-indexX;

         if(minX<x && x<maxX)
         {
           nnE[indexE][indexX]+=(1.0-coefE)*(1.0-coefX);
           nnE[indexE][indexX+1]+=(1.0-coefE)*coefX;
           nnE[indexE+1][indexX]+=coefE*(1.0-coefX);
           nnE[indexE+1][indexX+1]+=coefE*coefX;
         }
       }
       fclose(in);

       sprintf(str,"denE%d",step);
       out=fopen(str,"w");
       for(i=0; i<numE; i++)
       {
         for(j=0; j<numX; j++)
         {
           fprintf(out,"%d %g %g\n",i,j*dx+minX,nnE[i][j]*superP*1.6e-19);
         }
         fprintf(out,"\n");
       }
       fclose(out); 
       printf("denE%d is made.\n",step);

     }		//End of step
     for(i=0; i<=numE; i++)
       free(nnE[i]);
     free(nnE);
     break;

   case 2 :
     minX=atof(argv[3]);
     maxX=atof(argv[4]);
     numX=atoi(argv[5]);
     maxE=atof(argv[6]);
     step=atoi(argv[7]);
     superP=atof(argv[8]);

     dx=(maxX-minX)/((float)numX);

     numE=(int)(maxE);     

     nnE=(float **)malloc((numE+1)*sizeof(float *));
     for(i=0; i<=numE; i++)
       nnE[i]=(float *)malloc((numX+1)*sizeof(float ));

     for(i=0; i<numE; i++)
       for(j=0; j<numX; j++)
         nnE[i][j]=0.0;

     in = fopen(argv[2],"r");
//   fgets(str,100,in);
     while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&x,&y,&z,&ux,&uy,&uz,&gamma,&id,&core)!=EOF)
     {      
       energy=gamma*0.511;
       indexE=(int)(energy);
       coefE=energy-indexE;
       indexX=(int)((x-minX)/dx);
       coefX=(x-minX)/dx-indexX;

       if(minX<x && x<maxX)
       {
         nnE[indexE][indexX]+=(1.0-coefE)*(1.0-coefX);
         nnE[indexE][indexX+1]+=(1.0-coefE)*coefX;
         nnE[indexE+1][indexX]+=coefE*(1.0-coefX);
         nnE[indexE+1][indexX+1]+=coefE*coefX;
       }
     }
     fclose(in);

     sprintf(str,"denE%d",step);
     out=fopen(str,"w");
     for(i=0; i<numE; i++)
     {
       for(j=0; j<numX; j++)
       {
         fprintf(out,"%d %g %g\n",i,j*dx+minX,nnE[i][j]);
       }
       fprintf(out,"\n");
     }
     fclose(out); 
     printf("denE%d is made.\n",step);

     for(i=0; i<=numE; i++)
       free(nnE[i]);
     free(nnE);
     break;
   }
}
