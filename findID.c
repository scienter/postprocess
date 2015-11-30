#include "stdio.h"
#include "stdlib.h"
#include "math.h"
main(int argc, char *argv[])
{
   float x,y,z,ux,uy,uz,id,gamma,core;
   float minX,maxX,minP,maxP,kp,dt,v0,gamma0;
   FILE *in, *sample,*out;
   int mode,startCore,endCore,rank,i,cnt,dimension;
   int time,initial,final,timeStep,cores,laserFrame;
   char name[100];


   if(argc < 3)   { 
       printf("findID mode dimension\n");
       printf("mode 0 : file startCore endCore minX maxX minP maxP time laserFrame\n"); 
       printf("mode 1 : initial final timeStep cores laserFrame(using 'idSample' file)\n"); 
       exit(0);  
   }

   mode=atoi(argv[1]);
   dimension=atoi(argv[2]);

   kp=4.21e5;
   dt=1.67e-16;
   gamma0=14.931;
   v0=3e8*(1.0-0.5/gamma0/gamma0);

   if(mode==0)
   {
     startCore=atoi(argv[4]);
     endCore=atoi(argv[5]);
     minX=atof(argv[6]);
     maxX=atof(argv[7]);
     minP=atof(argv[8]);
     maxP=atof(argv[9]);
     time=atoi(argv[10]);
     laserFrame=atoi(argv[11]);

     switch (dimension)  {
     case 2:
       for(rank=startCore; rank<=endCore; rank++)
       {
         sprintf(name,"%s_%d",argv[3],rank);
         in = fopen(name,"r");
//     fgets(str,100,in);
         while(fscanf(in,"%g %g %g %g %g %g %g %g",&x,&y,&ux,&uy,&uz,&gamma,&id,&core)!=EOF)
         {
           if (x>minX && x<maxX && ux>minP  && ux<maxP)
           {
             if(laserFrame==0)
               printf("%g %g %g %g %g %g %.9g %g\n",x,y,ux,uy,uz,gamma,id,core);
             else if(laserFrame==1)
               printf("%g %g %g %g %g %g %.9g %g\n",x*kp-time*dt*v0*kp,y*kp,ux,uy,uz,gamma,id,core);
             else	;
           }
         }
         fclose(in);
       }
       break;
     case 3:
       for(rank=startCore; rank<=endCore; rank++)
       {
         sprintf(name,"%s_%d",argv[3],rank);
         in = fopen(name,"r");
//     fgets(str,100,in);
         while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&x,&y,&z,&ux,&uy,&uz,&gamma,&id,&core)!=EOF)
         {
           if (x>minX && x<maxX && ux>minP  && ux<maxP)
           {
             if(laserFrame==0)
               printf("%g %g %g %g %g %g %g %.9g %g\n",x,y,z,ux,uy,uz,gamma,id,core);
             else if(laserFrame==1)
               printf("%g %g %g %g %g %g %g %.9g %g\n",x*kp-time*dt*v0*kp,y*kp,z*kp,ux,uy,uz,gamma,id,core);
             else	;
           }
         }
         fclose(in);
       }
       break;
     }
   }

   else if(mode==1) 
   {
     initial=atoi(argv[3]);
     final=atoi(argv[4]);
     timeStep=atoi(argv[5]);
     cores=atoi(argv[6]);
     laserFrame=atoi(argv[7]);
   

     cnt=0;
     if(dimension==2)
     {
       sample=fopen("idSample","r");
       while(fscanf(sample,"%g %g %g %g %g %g %g %g %g"
                    ,&x,&y,&z,&ux,&uy,&uz,&gamma,&id,&core)!=EOF) 
         cnt=cnt+1;
       fclose(sample);
     }
     else if(dimension==3)
     {
       sample=fopen("idSample","r");
       while(fscanf(sample,"%g %g %g %g %g %g %g %g %g"
                    ,&x,&y,&z,&ux,&uy,&uz,&gamma,&id,&core)!=EOF) 
         cnt=cnt+1;
       fclose(sample);
     }
     float idList[cnt],coreList[cnt];

     switch (dimension) {
     case 2 :
       sample=fopen("idSample","r");
       for(i=0; i<cnt; i++)
         fscanf(sample,"%g %g %g %g %g %g %g %g"
                  ,&x,&y,&ux,&uy,&uz,&gamma,&idList[i],&coreList[i]); 
       fclose(sample);

       for(time=initial; time<=final; time+=timeStep)
       {
         sprintf(name,"id%d",time);
         out = fopen(name,"w");
         for(rank=0; rank<cores; rank++)
         {
           sprintf(name,"0Particle%d_%d",time,rank);
           in = fopen(name,"r");
           while(fscanf(in,"%g %g %g %g %g %g %g %g"
                      ,&x,&y,&ux,&uy,&uz,&gamma,&id,&core)!=EOF) 
           {
             for(i=0; i<cnt; i++) 
             {
               if (id==idList[i] && core==coreList[i]) 
                 fprintf(out,"%g %g %g %g %g %g %.9g %g\n",x,y,ux,uy,uz,gamma,id,core);
             }
           }
           fclose(in);
           printf("%s is done.\n",name);
         }	//End of rank  
         fclose(out);
       }		//End of time
       break;
     case 3 :
       sample=fopen("idSample","r");
       for(i=0; i<cnt; i++)
         fscanf(sample,"%g %g %g %g %g %g %g %g %g"
                  ,&x,&y,&z,&ux,&uy,&uz,&gamma,&idList[i],&coreList[i]); 
       fclose(sample);

       for(time=initial; time<=final; time+=timeStep)
       {
         if(laserFrame==0)
           sprintf(name,"id%d",time);
         else if(laserFrame==1)
           sprintf(name,"idLaser%d",time);
         else 	;
         out = fopen(name,"w");
         for(rank=0; rank<cores; rank++)
         {
           sprintf(name,"0Particle%d_%d",time,rank);
           in = fopen(name,"r");
           while(fscanf(in,"%g %g %g %g %g %g %g %g %g"
                      ,&x,&y,&z,&ux,&uy,&uz,&gamma,&id,&core)!=EOF) 
           {
             for(i=0; i<cnt; i++) 
             {
               if (id==idList[i] && core==coreList[i])
               { 
                 if(laserFrame==0)
                   fprintf(out,"%g %g %g %g %g %g %g %.9g %g\n",x,y,z,ux,uy,uz,gamma,id,core);
                 else if(laserFrame==1)
                   fprintf(out,"%g %g %g %g %g %g %g %.9g %g\n",x*kp-dt*time*v0*kp,y*kp,z*kp,ux,uy,uz,gamma,id,core);
                 else	;
               }
             }
           }
           fclose(in);
           printf("%s is done.\n",name);
         }	//End of rank  
         fclose(out);
       }		//End of time
       break;
     }
   }
}
