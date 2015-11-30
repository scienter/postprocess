#include "stdio.h"
#include "stdlib.h"
#include "math.h"
main(int argc, char *argv[])
{
   float x,y,z,ux,uy,uz,id,gamma,core,reCore,reId;
   float sumX,sumPx;
   FILE *in, *sample,*out;
   int step,cnt,dimension;
   int time,initial,final,timeStep,cores;
   char name[100];


   if(argc < 4)   { 
       printf("calPx dimension initial final timeStep \n");
       exit(0);  
   }

   dimension=atoi(argv[1]);
   initial=atoi(argv[2]);
   final=atoi(argv[3]);
   timeStep=atoi(argv[4]);

   if (dimension==2)
   {
     out = fopen("avePx","w");
     for(step=initial; step<=final; step+=timeStep)
     { 
       sprintf(name,"id%d",step);
       in = fopen(name,"r");
//     fgets(str,100,in);
       sumX=sumPx=0;
       cnt=0;
       while(fscanf(in,"%g %g %g %g %g %g %g %g",&x,&y,&ux,&uy,&uz,&gamma,&id,&core)!=EOF)
       {
         sumX+=x;
         sumPx+=ux;
         cnt++;
       }
       fclose(in);

       fprintf(out,"%d %g %g\n",step,sumX/cnt,sumPx/cnt);
       printf("id%d is done.\n",step);
     }
     fclose(out);
     printf("'avePx' file is made.\n");
   }
   else if (dimension==3)
   {
     out = fopen("avePx","w");
     for(step=initial; step<=final; step+=timeStep)
     { 
       sprintf(name,"id%d",step);
       in = fopen(name,"r");
//     fgets(str,100,in);
       sumX=sumPx=0;
       cnt=0;
       while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&x,&y,&z,&ux,&uy,&uz,&gamma,&id,&core)!=EOF)
       {
         sumX+=x;
         sumPx+=ux;
         cnt++;
       }
       fclose(in);

       fprintf(out,"%d %g %g\n",step,sumX/cnt,sumPx/cnt);
       printf("id%d is done.\n",step);
     }
     fclose(out);
     printf("'avePx' file is made.\n");
   }
}
