#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
void main(int argc, char *argv[])
{
  float minX,maxX,minE,maxE,dx,eenergy;
  float xx,tmp,gamma;
  int i,j,numParticle,numX,numE,time,indexI,indexJ;
  int cnt,min,max,remain,numData,rnk,cores;
  float **n,*E,*x;
  FILE *in,*out;
  char fileName[100];
  void restoreData();

  if(argc < 7)   { 
    printf("energy3D cores minX maxX numX minE maxE time\n");
    exit(0);  
  }

  cores=atoi(argv[1]);
  minX=atof(argv[2]);
  maxX=atof(argv[3]);
  numX=atoi(argv[4]);
  minE=atof(argv[5]);
  maxE=atof(argv[6]);
  time=atoi(argv[7]);
  

  dx=(maxX-minX)/((float) numX);
  numE=(int)((maxE-minE)/1e6);
  x=(float *)malloc(numX*sizeof(float ));
  E=(float *)malloc(numE*sizeof(float ));
  n=(float **)malloc(numX*sizeof(float *));
  for(i=0; i<numX; i++)
    n[i]=(float *)malloc(numE*sizeof(float ));
  for(i=0; i<numX; i++)
    for(j=0; j<numE; j++)
      n[i][j]=0.0;

  for(rnk=0; rnk<cores; rnk++)
  {
    sprintf(fileName,"0Particle%d_%d",time,rnk);
    in=fopen(fileName,"r");
    while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&xx,&tmp,&tmp,&tmp,&tmp,&tmp,&gamma,&tmp,&tmp)!=EOF)
    {
      indexI=(int)((xx-minX)/dx);
      eenergy=(gamma-1.0)*0.511*1e6;
      indexJ=(int)((eenergy-minE)/1e6);
      if(indexI>=0 && indexI<numX && indexJ>=0 && indexJ<numE)
      {
        n[indexI][indexJ]+=1.0;
      }
    }
    fclose(in);
    printf("rank%d is done.\n",rnk);
  }

  for(i=0; i<numX; i++)
    x[i]=i*dx+minX;
  for(i=0; i<numE; i++)
    E[i]=i*1e6+minE;

  sprintf(fileName,"denE%d",time);
  out=fopen(fileName,"w");
    for(i=0; i<numX; i++)
    {
      for(j=0; j<numE; j++)
        fprintf(out,"%g %g %g\n",x[i],E[j],n[i][j]);
      fprintf(out,"\n");
    }

  free(x);
  free(E);
  for(i=0; i<numX; i++)
    free(n[i]);
  free(n);

}
