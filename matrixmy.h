#include <stdio.h>
#include <stdlib.h>

#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
   fprintf(stderr,"Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   exit(1);
};
/*-----------------------------------------------------------------*/
double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{ double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in dvector()"); return v-nl+NR_END;
}
/*---------------------------------------------------------------------*/
double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{ long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
};
/*-----------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
long **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a long matrix with subscript range m[nrl..nrh][ncl..nch] */
{ long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  long **m;

  /* allocate pointers to rows */
  m=(long **) malloc((size_t)((nrow+NR_END)*sizeof(long*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(long *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(long)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
};
/*-----------------------------------------------------------------*/
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{ long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;

  /* allocate pointers to pointers to rows */
  t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
  if (!t) nrerror("allocation failure 1 in d3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
  if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");

  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
     t[i]=t[i-1]+ncol;
     t[i][ncl]=t[i-1][ncl]+ncol*ndep;
     for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
};
/*-------------------------------------------------------------------------*/
/*-------------------------my----------------------------------------*/
long ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a int 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{ long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  long ***t;

  /* allocate pointers to pointers to rows */
  t=(long ***) malloc((size_t)((nrow+NR_END)*sizeof(long**)));
  if (!t) nrerror("allocation failure 1 in i3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(long **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(long*)));
  if (!t[nrl]) nrerror("allocation failure 2 in i3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(long *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(long)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in i3tensor()");

  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
     t[i]=t[i-1]+ncol;
     t[i][ncl]=t[i-1][ncl]+ncol*ndep;
     for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
};
/*-------------------------------------------------------------------------*/
void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
   free((FREE_ARG) (v+nl-NR_END));
}
/*-----------------------------------------------------------------------*/
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}
/*-----------------------------------------------------------------------*/
void free_imatrix(long **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by imatrix() */
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}
/*-----------------------------------------------------------------------*/
void free_d3tensor(double ***t, long nrl, long nrh,
                   long ncl, long nch, long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
   free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
   free((FREE_ARG) (t[nrl]+ncl-NR_END));
   free((FREE_ARG) (t+nrl-NR_END));
}
/*------------------------my-----------------------------------------------*/
void free_i3tensor(long ***t, long nrl, long nrh,
                   long ncl, long nch, long ndl, long ndh)
/* free a long i3tensor allocated by i3tensor() */
{
   free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
   free((FREE_ARG) (t[nrl]+ncl-NR_END));
   free((FREE_ARG) (t+nrl-NR_END));
}
