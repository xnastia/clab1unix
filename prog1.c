#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define frand() ((double) rand() / (RAND_MAX));
#include "matrixmy.h"

//#define N 3
//#define It 2
#define eps 1e-8
FILE *fp;
FILE *fp_ini;
FILE *fp_x;
FILE *fp_v;
FILE *fp_orient;
FILE *fp_vcm;
FILE *fp_pcm;
FILE *fp_b;


float dist(float x1, float x2, float y1, float y2){
  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

int main(void)
/* 0 */{int i,j,k,n;
long N; long It;
float R,dt,d,delta_max,delta,delta0,Tcalc;

//	   printf("\nPlease, input the average number of neighbors \n");
//	   printf("n = ");
//	   scanf("%i", &n);
 
//	   printf("\n Please, input the number of particles N\n");
//       printf("N = ");
//       scanf("%li", &N);

//	   printf("\n Please, input the number of iterations It\n");
//       printf("It = ");
//       scanf("%li", &It);

/*	   printf("\n Only every m-th step will be recorded\n");
       printf("m = ");
       scanf("%i", &m);*/

//	   printf("\nPlease, input the noise amplitude delta_max (in degrees)\n");
//       printf("delta_max = ");
//       scanf("%f", &delta_max);
 n = 100;
 N = 10000;
 It = 1000;
 delta_max = 10;
	   delta0 = delta_max/180.*3.14159;
	   R =pow(n/(N*3.14159),0.5);
	   dt = 2./It;

printf("\n\nOK. The calculation will be performed ");
printf("for N = %li particles \n in It = %li iterations \n", N,It);
//printf("with every m = %i -th step recorded \n", m);
printf("with time increment dt=%f , interaction range R=%f\n",dt,R);
printf("and the average number within the interaction range n=%i\n",n);
printf("and the noise opening range delta=%f\n",2*delta_max);

		fp = fopen("dynamic2dim.txt", "w");
		fp_ini = fopen("ini.txt", "w");
		fp_vcm = fopen("vcm2dim.txt", "w");
		fp_pcm = fopen("pcm2dim.txt", "w");
		fp_x = fopen("x2dim.txt", "w");
		fp_v = fopen("v2dim.txt", "w");
		fp_orient = fopen("orient2dim.txt","w");

fprintf(fp,"The calculation is performed with the following parameters\n");
fprintf(fp,"N  = %li\n",N);
fprintf(fp,"It = %li\n",It);
//fprintf(fp,"m = %i\n",m);
fprintf(fp,"R  = %f\n",R);
fprintf(fp,"dt = %f\n",dt);
fprintf(fp,"neighb = %i\n",n);
fprintf(fp,"delta_max = %f\n",delta_max);
                                                                                                                                                                                                                                                                                                                                             
fprintf(fp_ini,"%li\n",N);
fprintf(fp_ini,"%li\n",It);
//fprintf(fp_ini,"%i\n",m);
fprintf(fp_ini,"%f\n",R);
fprintf(fp_ini,"%f\n",dt);
fprintf(fp_ini,"%i\n",n);
fprintf(fp_ini,"%f\n",delta_max);

/* here should start the external loop for the initial config */
	double sx,sy,sx1,sy1,s_mod,ro,rho,ca,sa,vx1,vy1,r,dt1,x1,y1;
	double*   angle; /* angle[]     */
	//double**   orient; /* orientation[][] */
	double**  x; /* x[][]   */
	double**  y; /* y[][][] */
	double*  xc; /* x[][]   */
	double*  yc; /* y[][][] */
 	double**  vx;/* vx[][][] */
 	double**  vy;/* vx[][][] */
 	double*  vx_centerofmass;
 	double*  vy_centerofmass;
 	double*  vmod_centerofmass;
//	double*  v0;
// 	long**  B; /* B[][] */

 angle=dvector(0,N-1);
 vx_centerofmass = dvector(0,It);
 vy_centerofmass = dvector(0,It);
 xc = dvector(0,It);
 yc = dvector(0,It);
 vmod_centerofmass = dvector(0,It);
 x=dmatrix(0,N-1,0,It);
 y=dmatrix(0,N-1,0,It);
 vx=dmatrix(0,N-1,0,It);
 vy=dmatrix(0,N-1,0,It);
//	v0=dvector(0,N);
// B = imatrix(0,N,0,N);

//for(n=0;n<=N_cfg;n++){ //This is the generating of the initial config

srand((unsigned)time(NULL));
sx=0;sy=0;
for( i=0;i<=N-1;i++)
/* 1 */ {
ro =frand();ro=sqrt(ro);
// the initial (random) coordinates and velocities of the boids

/*in a box*/
sx=1-2*frand();
sy=1-2*frand();//
x[i][0]=sx;y[i][0]=sy;

/*v0[i] = 0.1+0.9*frand();
angle[i] = 1.-2.*frand();
x[i][0]=ro*cos(angle[i]*asin(1)*2.);
y[i][0]=ro*sin(angle[i]*asin(1)*2.);*/

// the initial ranfom direction of the boid's velocity
angle[i] = 1.-2.*frand();

vy[i][0]=sin(angle[i]*asin(1)*2.); // the initial velocity components
vx[i][0]=cos(angle[i]*asin(1)*2.);
//printf("initial coordinates and velocity of boid No. %i\n",i+1);

fprintf(fp_x,"% .4f\t% .4f\n",x[i][0],y[i][0]);
fprintf(fp_v,"% .4f\t% .4f\n",vx[i][0],vy[i][0]);

/* 1 */ }


/*---------------Center of mass velocity--------------------*/
sx=0;sy=0;
for(i=0;i<=N-1;i++)
{
sx+= vx[i][0];
sy+= vy[i][0];
}

vx_centerofmass[0] = sx/N; vy_centerofmass[0] = sy/N;
vmod_centerofmass[0] = sqrt(pow(vx_centerofmass[0],2)+pow(vy_centerofmass[0],2));
fprintf(fp_vcm,"% .6f\t% .6f\t% .6f\n",
vx_centerofmass[0],vy_centerofmass[0],vmod_centerofmass[0]);

sx=0;sy=0;
for(i=0;i<=N-1;i++)
{
sx+= x[i][0];
sy+= y[i][0];
}

xc[0] = sx/N; yc[0] = sy/N;
fprintf(fp_pcm,"% .6f\t% .6f\n", xc[0],yc[0]);

/*---------------------------------------------------------------*/
free_dvector(angle,0,N-1);

/*----------------------------Main body--------------------------*/

printf("\n\nHere comes the evolution\n\n");
/*----------------------------Main body------------------------*/
 

for( k = 1;k<=It;k++) //start iteration loop
/* 2 */{
for(i=0;i<=N-1;i++)		//start boid loop
/* 3 */{sx=0; sy=0;

sx=x[i][k-1]+vx[i][k-1]*dt;
sy=y[i][k-1]+vy[i][k-1]*dt;

/* updating the position with cyclic boundary conditions */
if(fabs(sx)>1.) sx1=1.; else sx1=0;
if(fabs(sy)>1.) sy1=1.; else sy1=0;
x[i][k]= sx-2.*sx/fabs(sx)*sx1;
y[i][k]= sy-2.*sy/fabs(sy)*sy1;

sx=0; sy=0;

/*----------Calculating the distance and incidence matrices-------*/
//for( j=0;j<=N-1;j++)
///* 4 */ {						//start of boid's neighbor loop
//if(i!=j &&
//dist(x[i][k-1],x[j][k-1],y[i][k-1],y[j][k-1]) <= R) B[i][j]= 1;
//else B[i][j] = 0;
///* 4 */}						// end of boid's neighbor loop
/*---------------------------------------------------------------*/


/*------------------we continue with the velocity---------------------*/
for( j=0;j<=N-1;j++)
/* 5 */{						//start of boid's velocity loop
if(i!=j && dist(x[i][k-1],x[j][k-1],y[i][k-1],y[j][k-1]) <= R)
{sx+= (vx[j][k-1]);
sy+= (vy[j][k-1]);}
else 
				{sx+= 0; sy+= 0;}
/* 5 */ 
	   }

if(delta0 == 0){
//sx=sx;sy=sy;						//end of boid's velocity loop
s_mod = sqrt(pow(sx,2)+pow(sy,2));
if(s_mod!=0) {vx[i][k]=sx/s_mod;vy[i][k]=sy/s_mod;}
else {vx[i][k]=vx[i][k-1];vy[i][k]=vy[i][k-1];}
}
			else
//noise
delta = 1.-2.*frand();
sx=sx*cos(delta*delta0)-sy*sin(delta*delta0);			
sy=sy*cos(delta*delta0)+sx*sin(delta*delta0);
//
s_mod = sqrt(pow(sx,2.)+pow(sy,2.));
if(s_mod!=0) {vx[i][k]=sx/s_mod;vy[i][k]=sy/s_mod;}
else {vx[i][k]=vx[i][k-1]*cos(delta*delta0)-vy[i][k-1]*sin(delta*delta0);
				vy[i][k]=vy[i][k-1]*cos(delta*delta0)+vx[i][k-1]*sin(delta*delta0);}


//if(k%m==0){
fprintf(fp_x,"% .4f\t% .4f\n",x[i][k],y[i][k]);
fprintf(fp_v,"% .4f\t% .4f\n",vx[i][k],vy[i][k]);
//}

////////////writing data to file////////////////////////////////
						// end of boid loop
/*3*/}
/*---------------Center of mass velocity--------------------*/
sx=0;sy=0;
for(i=0;i<=N-1;i++)
{
sx+= vx[i][k];
sy+= vy[i][k];
}

vx_centerofmass[k] = sx/N; vy_centerofmass[k] = sy/N;
vmod_centerofmass[k] = sqrt(pow(vx_centerofmass[k],2)+pow(vy_centerofmass[k],2));
fprintf(fp_vcm,"% .6f\t% .6f\t% .6f\n",
vx_centerofmass[k],vy_centerofmass[k],vmod_centerofmass[k]);
/*---------------------------------------------------------*/

/*---------------Center of mass position--------------------*/
sx=0;sy=0;
for(i=0;i<=N-1;i++)
{
sx+= x[i][k];
sy+= y[i][k];
}

xc[k] = sx/N; yc[k] = sy/N;
fprintf(fp_pcm,"% .6f\t% .6f\n",
xc[k],yc[k]);
/*---------------------------------------------------------*/

//fprintf(fp,"\n");

//if(fabs(vmod_centerofmass[k]- vmod_centerofmass[k-1]) <= eps && //
//fabs(vmod_centerofmass[k]-vmod_centerofmass[k-1])>=1e-10)//
//{fprintf(fp,"The fixed configuration is reached at %i-th step t=%f\n",k,k*dt); goto stp;} else;

/* 2 */}						// end of iteration loop

/////////////////////////////////////////
	//making orientation array 
	for( k=0;k<=It;k++)
		{
		for (i=0;i<=N-1;i++)
			{sx=0;
				sx = (1-(vx[i][k]*vx[i][It]+vy[i][k]*vy[i][It]) ) / (1+(vx[i][k]*vx[i][It]+vy[i][k]*vy[i][It]) );
	fprintf(fp_orient,"%.4f\t",sx);		
		}//end of boid loop
	fprintf(fp_orient,"\n");		
	}//end of iteration loop
////////////////////////////////////////
	goto stp;// comment this line to reach the code below

stp:

 free_dvector(vx_centerofmass,0,N);
 free_dvector(vy_centerofmass,0,N);
 free_dvector(vmod_centerofmass,0,N);
 free_dmatrix(x,0,N,0,It);
 free_dmatrix(y,0,N,0,It);
 free_dvector(xc,0,N);
 free_dvector(yc,0,N);
 free_dmatrix(vx,0,N,0,It);
 free_dmatrix(vy,0,N,0,It);
 //free_imatrix(B,0,N,0,N);

 //}
Tcalc=clock()/CLOCKS_PER_SEC;

	printf ("\n It took %.1f seconds.\n", Tcalc);
	printf("\n End of program \n");

   return 0;
/* 0 */ }

