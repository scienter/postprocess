// From x-y ascii data from xgrafix, reconstruct the density profile
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

void main(int argc, char *argv[])
{
   float R,centerX,centerY,x,y;
   int i;

   if(argc < 3)
   {
      printf("makecircle R centerX centerY\n");
      exit(0);
   }
   R=atof(argv[1]);
   centerX=atof(argv[2]);
   centerY=atof(argv[3]);

   for (i=0; i<360; i++)
   {
     x=centerX+R*cos(i*3.14/180.0);
     y=centerY+R*sin(i*3.14/180.0);
     printf("%g %g %g\n",x,y,0.0);
   }
}
