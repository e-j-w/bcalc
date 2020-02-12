#include "bcalc.h"

void printHelp(){
  printf("B(E2) Calculator -- calculates reduced transition probabilities.\n");
  printf("example usage: bcalc -e VALUE -m VALUE -lt VALUE\n");
  printf("  Both of the following are needed:\n");
  printf("    -e         --  transition energy in keV\n");
  printf("    -m         --  multipole (eg. E1, M1, E2, etc.)\n");
  printf("  One of the following is needed:\n");
  printf("    -lt        --  Mean transition lifetime (in ps)\n");
  printf("    -hl        --  Transition half-life (in ps)\n");
  printf("    -b         --  Reduced transition probability (for the\n"); 
  printf("                   L multipole) in units of e^2 fm^(2L) for\n");
  printf("                   electric multipoles or uN^2 fm^(2L-2) for\n");
  printf("                   magnetic multipoles.\n");
  printf("  Optional parameters:\n");
  printf("    -d         --  mixing ratio with the L+1 multipole\n");
  printf("    -ji        --  inital spin (integer or half-integer)\n");
  printf("    -jf        --  final spin (integer or half-integer)\n");
  printf("  Flags:\n");
  printf("    --up       --  Use/calculate transition probability from\n");
  printf("                   final to initial state instead of vice versa.\n");
  printf("                   If used, requires the -ji and -jf parameters.\n");
  printf("    --verbose  --  Print detailed information.\n");
}

double dblfac(unsigned int n){ 
  double val = 1.0;
  int i;
  for (i=n; i>=0; i=i-2){
    if (i==0 || i==1)
      return val;
    else
      val *= i;
  }
  return val;
}

/* calcualtes the value of the reduced transtion probability */
void calcB(int bup, int EM, int L, double Et, double lt, double ji, double jf, int verbose, char * mstr){
  if(verbose){
    printf("\nB(%s) CALCULATION\n-----------------\n",mstr);
  }

  double fac = 8.0*PI*(L+1)/(L*HBAR_MEVS*pow(dblfac((2.0*L)+1.0),2.0)*LN2);
  fac = fac * (pow((Et/HBARC_MEVFM),2.0*L + 1.0));
  /* lifetime to reduced transition probability */
  if(EM==0){
    /* electric */
    /* e^2 fm^(2L) */
    double b=1/(fac*lt);
    if(bup){
      b = b*(2.0*ji + 1.0)/(2.0*jf + 1.0);
    }
    if(L>0){
      printf("%0.4E e^2 fm^%i\n",b,2*L);
      printf("%0.4E e^2 b^%i\n",b/pow(BARN_FM,L),L);
    }else{
      printf("%0.4E e^2\n",b);
    }
    
  }else if(EM==1){
    /* magnetic */
    /* uN^2 fm^(2L-2) */
    fac = fac * UN_MEVFM3/ESQ_MEVFM;
    //printf("fac: %f\n",fac);
    double b=1/(fac*lt);
    if(bup){
      b = b*(2.0*ji + 1.0)/(2.0*jf + 1.0);
    }
    if(L>1)
      printf("%0.4E uN^2 fm^%i\n",b,(2*L) - 2);
    else
      printf("%0.4E uN^2\n",b);
  }
}

void calcLt(int bup, int EM, int L, double Et, double b, double ji, double jf, int verbose, char * mstr){
  if(verbose){
    printf("\nLIFETIME CALCULATION\n--------------------\n");
  }

  if(bup){
    b = b*(2.0*jf + 1.0)/(2.0*ji + 1.0);
  }

  double fac = 8.0*PI*(L+1)/(L*HBAR_MEVS*pow(dblfac((2.0*L)+1.0),2.0)*LN2);
  fac = fac * (pow((Et/HBARC_MEVFM),2.0*L + 1.0));
  /* reduced transition probability to lifetime */
  if(EM==1){
    fac = fac * UN_MEVFM3/ESQ_MEVFM;
  }
  double lt = 1/(fac*b);
  lt = lt / 1.0E-12; //convert lifetime to ps
  printf("%0.4E ps\n",lt);
}

int main(int argc, char *argv[]) {

  if (argc == 1) {
    printHelp();
    exit(-1);
  }

  int i; /*counters*/

  /*initialize parameter values*/
  char mstr[2], mstr1[2];
  double Et = -1.; /* transition energy */
  int L = -1; /* multipolarity */
  int EM = -1; /* 0=electric, 1=magnetic */
  double lt = 0.; /* lifetime */
  double lt1 = 0.; /* lifetime for L+1 multipole (mixed transitions) */
  double b = 0.; /* reduced transition probability */
  int calcMode = -1; /* 0=calulate B, 1=calculate lifetime */
  int verbose = 0; /* 0=none, 1=verbose */
  int bup = 0; /* 0=down, 1=up */
  double ji = -1.; /* Initial spin */
  double jf = -1.; /* Final spin*/
  double delta = 0; /* mixing ratio */
  int useDelta = 0; /* 0=no mixing, 1=mixing */

  /*read parameters*/
  for(i=0;i<argc;i++){
    if(strcmp(argv[i],"--verbose")==0){
      verbose = 1;
    }else if(strcmp(argv[i],"--up")==0){
      bup = 1;
    }else if(strcmp(argv[i],"--help")==0){
      printHelp();
      exit(-1);
    }
  }
  for(i=0;i<(argc-1);i++){
    if((strcmp(argv[i],"-E")==0)||(strcmp(argv[i],"-e")==0)){
      Et=atof(argv[i+1]);
      if(Et <= 0.){
        printf("ERROR: Invalid transition energy (%f).  Value must be a positive number.\n",Et);
        exit(-1);
      }
    }else if((strcmp(argv[i],"-M")==0)||(strcmp(argv[i],"-m")==0)){
      strncpy(mstr,argv[i+1],sizeof(mstr));
      if(mstr[0] == 'E'){
        EM=0;
      }else if(mstr[0] == 'M'){
        EM=1;
      }else{
        printf("ERROR: invalid multipole value.\n");
        exit(-1);
      }
      if(isdigit(mstr[1])!=0){
        L=mstr[1] - '0';
      }else{
        printf("ERROR: invalid multipole value.\n");
        exit(-1);
      }
      
    }else if((strcmp(argv[i],"-Lt")==0)||(strcmp(argv[i],"-lt")==0)){
      lt=atof(argv[i+1]);
      calcMode = 0;
    }else if((strcmp(argv[i],"-Hl")==0)||(strcmp(argv[i],"-hl")==0)){
      lt=atof(argv[i+1])/LN2;
      calcMode = 0;
    }else if((strcmp(argv[i],"-B")==0)||(strcmp(argv[i],"-b")==0)){
      b=atof(argv[i+1]);
      calcMode = 1;
    }else if(strcmp(argv[i],"-d")==0){
      delta=atof(argv[i+1]);
      useDelta = 1;
    }else if(strcmp(argv[i],"-ji")==0){
      ji=atof(argv[i+1]);
      if(ji<0){
        printf("ERROR: Invalid initial spin value provided (%0.1f).  The value must be a positive integer or half-integer.\n",ji);
        exit(-1);
      }
    }else if(strcmp(argv[i],"-jf")==0){
      jf=atof(argv[i+1]);
      if(jf<0){
        printf("ERROR: Invalid final spin value provided (%0.1f).  The value must be a positive integer or half-integer.\n",jf);
        exit(-1);
      }
    }
  }

  /*check argument values for validity*/
  if(Et < 0.){
    printf("ERROR: Missing parameter.\n");
    printf("    -e  --  transition energy in keV\n");
    exit(-1);
  }
  if(calcMode < 0){
    printf("ERROR: One of the following parameters is needed:\n");
    printf("   -lt  --  Mean transition lifetime (in ps)\n");
    printf("   -hl  --  Transition half-life (in ps)\n");
    printf("   -b   --  Reduced transition probability\n");
    exit(-1);
  }
  if((EM==1) && (L ==0)){
    printf("Magnetic monopole transitions are not allowed.\n");
    exit(-1);
  }
  if(L<0){
    printf("ERROR: Missing parameter.\n");
    printf("    -m  --  multipole (eg. E1, M1, E2, etc.)\n");
    exit(-1);
  }
  if ((fmod(ji,0.5)!=0.) || (fmod(jf,0.5)!=0.) || (fmod(jf-ji,1)!=0.)){
     printf("ERROR: Initial and final spins must both be either integer or half-integer.\n");
     exit(-1);
  }
  if((bup == 1)&&((ji == -1)||(jf == -1))){
    printf("WARNING: To calculate B(%s) up, the initial and final spin must be known.\nAssuming a 2 -> 0 transition.\n",mstr);
    ji=2.;
    jf=0.;
  }

  /*print extra info*/
  if(verbose){
    printf("\nINPUT PARAMETERS\n----------------\n");
    printf("Transition energy: %0.3f keV\n",Et);
    printf("Transition multipole: ");
    if(EM==0)
      printf("electric ");
    else
      printf("magnetic ");
    if(L==0)
      printf("monopole");
    else if(L==1)
      printf("dipole");
    else if(L==2)
      printf("quadrupole");
    else if(L==3)
      printf("octopole");
    else if(L==4)
      printf("hexadecapole");
    else
      printf("L = %i",L);
    if(useDelta&&(calcMode==0)){
      printf(" (L+1 mixing, delta = %0.3f)", delta);
    }
    printf("\n");
    if(calcMode == 0)
      printf("Mean lifetime: %0.3f ps\n",lt);
    else if(calcMode == 1){
      printf("B(%s): %f ",mstr, b);
      if(EM==0){
        printf("e^2 fm^%i\n",2*L);
      }else if(EM==1){
        printf("uN^2 fm^%i\n",(2*L) - 2);
      }
    }
  }

  /* mixing ratio calculation */
  if(useDelta){
    lt1 = lt * (1.0 + delta*delta) / (delta*delta);
    lt = lt * (1.0 + delta*delta);
    if(mstr[0] == 'E')
      sprintf(mstr1,"M%i",L+1);
    else
      sprintf(mstr1,"E%i",L+1);
    
  }

  lt=lt*1.0E-12; //convert lifetime to s
  Et=Et/1000.0; //convert energy to MeV

  if(calcMode == 0){
    calcB(bup,EM,L,Et,lt,ji,jf,verbose,mstr);
    if(useDelta)
      calcB(bup,!EM,L+1,Et,lt1,ji,jf,verbose,mstr1);
  }else if(calcMode == 1){
    calcLt(bup,EM,L,Et,b,ji,jf,verbose,mstr);
  }
  

  return 0;
}