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

char mstr[2];
double Et; /* transition energy */
int L; /* multipolarity */
int EM; /* 0=electric, 1=magnetic */
double lt; /* lifetime */
double b; /* reduced transition probability */
int calcMode; /* 0=calulate B, 1=calculate lifetime */
int verbose; /* 0=none, 1=verbose */
int bup; /* 0=down, 1=up */
int ji,jf; /* Initial and final spins*/