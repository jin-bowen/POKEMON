#define TRUE  1
#define FALSE 0
typedef int BOOL;
#define UseDouble 1             /* all floating point double */

#ifdef UseDouble
   typedef double real;
#else
   typedef float real;
#endif


real qf_swig(real*,int,real*,int,int*,int,real,real,int,real,real*,int,int*,int);

real qf(real*,real*,int*,int,real,real,int,real,real*,int*);
