#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define pickY 1

int main(int argc, char *argv[])
{
   float x,y,z,Ex,Ey,Ez,Bx,By,Bz, pickValY,pickValZ,tmp;
   float dev;
   FILE *in;
   char str[100];
   int whatMode();
   int i,mode,dimension;

   if(argc < 2)   { 
      printf("pickField mode dimension\n");
      printf("mode -> pickY File pickValY pickValZ dev\n");
      exit(0);
   }
   mode = whatMode(argv[1]);
   dimension = atoi(argv[2]);

   switch((mode-1)*3+dimension)   {
   case (pickY-1)*3+2 :
     pickValY = atof(argv[4]);
     pickValZ = atof(argv[5]);
     dev= atof(argv[6]);

     in = fopen(argv[3],"r");
//      fgets(str,100,inCen);       
     while(fscanf(in,"%g %g %g %g %g %g %g %g",&x,&y,&Ex,&Ey,&Ez,&Bx,&By,&Bz)!=EOF) 
     {
       if (y>=pickValY-dev && y<pickValY+dev)         
         printf("%g %g %g %g %g %g %g\n",x,Ex,Ey,Ez,Bx,By,Bz);
     }
     fclose(in);
     break;
   case (pickY-1)*3+3 :
     pickValY = atof(argv[4]);
     pickValZ = atof(argv[5]);
     dev= atof(argv[6]);

     in = fopen(argv[3],"r");
//      fgets(str,100,inCen);       
     while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&x,&y,&z,&Ex,&Ey,&Ez,&Bx,&By,&Bz)!=EOF) 
     {
       if (y>=pickValY-dev && y<pickValY+dev &&         
           z>=pickValZ-dev && z<pickValZ+dev)          
         printf("%g %g %g %g %g %g %g\n",x,Ex,Ey,Ez,Bx,By,Bz);
     }
     fclose(in);
     break;
   default :
     printf("I could not find a proper mode!.\n");
   }
}

int whatMode(char *str)
{
  if(strstr(str, "pickY"))	return pickY;
  else				return 1e5;
}
