#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define PI         3.14159265359
#define LN2        0.693147180559945
#define ESQ_MEVFM      1.44 /* MeV fm */
#define UN_MEVFM3      0.015922 /* MeV fm^3 */
#define UNSQ           1.5922E-38 /* keV cm^3 */
#define BARN_FM        100.0 /* fm^2 */
#define HBAR_MEVS      6.58212E-22 /* MeV s */
#define HBARC_MEVFM    197.327 /* MeV fm */

/* function prototypes */
void printHelp(void);
double dblfac(unsigned int);
double ltsp(const int,const int,const int,const double);
void calcBeta2(const double,const double,const int,const int,const int,const int);
void calcBeta2Lt(const double,const double,const int,const int,const int);
void calcB(const int,const int,const int,const double,const double,const double,const double,const int,const int,const char *,const int);
void calcLt(const int,const int,const int,const double,double,const double,const double,const int,const int,const char *,const int,const double);
