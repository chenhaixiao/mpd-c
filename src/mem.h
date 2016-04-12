#ifndef __mem_h__
#define __mem_h__
#endif

#ifndef UINT_TYPE
typedef unsigned int uint;
typedef unsigned char uchar;
#define UINT_TYPE 1
#endif

int       *ivector(int, int);
char      *cvector(int, int);
uchar     *ucvector(int, int);
double    *dvector(int, int);
double   **dmatrix(int, int, int, int);
char     **cmatrix(int, int, int, int);
uchar    **ucmatrix(int, int, int, int);
void       free_cvector(char *, int, int);
void       free_ucvector(uchar *, int, int);
void       free_ivector(int *, int, int);
void       free_dvector(double *, int, int);
void       free_dmatrix(double **, int, int, int, int);
void       free_cmatrix(char **, int, int, int, int);
void       free_ucmatrix(uchar **, int, int, int, int);
int      **imatrix(int, int, int, int);
void       free_imatrix(int **, int, int, int, int);
