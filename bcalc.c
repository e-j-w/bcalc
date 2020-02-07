#include "bcalc.h"

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

int main(int argc, char *argv[]) {

  if (argc == 1) {
    printf("B(E2) Calculator -- calculates reduced transition probabilities.\n");
    printf("usage: bcalc -E VALUE -M VALUE -Lt VALUE -Hl VALUE -B VALUE\n");
    printf("   Both of the following are needed:\n");
    printf("      -E         --  transition energy in keV\n");
    printf("      -M         --  multipole (eg. E1, M1, E2, etc.)\n");
    printf("   One of the following is needed:\n");
    printf("      -Lt        --  Mean transition lifetime (in ps)\n");
    printf("      -Hl        --  Transition half-life (in ps)\n");
    printf("      -B         --  Reduced transition probability\n");
    printf("   Optional parameters:\n");
    printf("      -up        --  Use/calculate transition probability from\n");
    printf("                     final to initial state instead of vice versa.\n");
    printf("                     If used, requires the parameters:\n");
    printf("                         -ji   --  inital spin\n");
    printf("                         -jf   --  final spin\n");
    printf("      -verbose   --  Print extra info\n");
    exit(-1);
  }

  int i; /*counters*/

  /*initialize parameter values*/
  Et = -1.;
  L = -1;
  lt = 0.;
  b = 0.;
  calcMode = -1;
  verbose = 0;
  bup = 0;
  ji = -1;
  jf = -1;

  /*read parameters*/
  for(i=0;i<argc;i++){
    if(strcmp(argv[i],"-verbose")==0){
      verbose = 1;
    }else if(strcmp(argv[i],"-up")==0){
      bup = 1;
    }
  }
  for(i=0;i<(argc-1);i++){
    if(strcmp(argv[i],"-E")==0){
      Et=atof(argv[i+1]);
      if(Et <= 0.){
        printf("ERROR: Invalid transition energy (%f).  Value must be a positive number.\n",Et);
        exit(-1);
      }
    }else if(strcmp(argv[i],"-M")==0){
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
      
    }else if(strcmp(argv[i],"-Lt")==0){
      lt=atof(argv[i+1]);
      calcMode = 0;
    }else if(strcmp(argv[i],"-Hl")==0){
      lt=atof(argv[i+1])/LN2;
      calcMode = 0;
    }else if(strcmp(argv[i],"-B")==0){
      b=atof(argv[i+1]);
      calcMode = 1;
    }else if(strcmp(argv[i],"-ji")==0){
      ji=atoi(argv[i+1]);
      if(ji<0){
        printf("ERROR: Invalid initial spin value provided (%i).  The value must be a positive integer.\n",ji);
        exit(-1);
      }
    }else if(strcmp(argv[i],"-jf")==0){
      jf=atoi(argv[i+1]);
      if(jf<0){
        printf("ERROR: Invalid final spin value provided (%i).  The value must be a positive integer.\n",jf);
        exit(-1);
      }
    }
  }

  /*check argument values*/
  if(Et < 0.){
    printf("ERROR: Missing parameter.\n");
    printf("    -E  --  transition energy in keV\n");
    exit(-1);
  }
  if(calcMode < 0){
    printf("ERROR: One of the following parameters is needed:\n");
    printf("   -Lt  --  Mean transition lifetime (in ps)\n");
    printf("   -Hl  --  Transition half-life (in ps)\n");
    printf("   -B   --  Reduced transition probability\n");
    exit(-1);
  }
  if((EM==1) && (L ==0)){
    printf("Magnetic monopole transitions are not allowed.\n");
    exit(-1);
  }
  if((bup == 1)&&((ji == -1)||(jf == -1))){
    printf("ERROR: To calculate B(%s) up, the initial and final spin must be specified.\n",mstr);
    exit(-1);
  }

  /*print extra info*/
  if(verbose){
    printf("INPUT PARAMETERS\n----------------\n");
    printf("Transition energy: %f keV\n",Et);
    printf("Transition multipolarity: ");
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
    printf("\n");
    if(calcMode == 0)
      printf("Mean lifetime: %f ps\n",lt);
    else if(calcMode == 1){
      printf("B(%s): %f ",mstr, b);
      if(EM==0){
        printf("e^2 fm^%i\n",2*L);
      }else if(EM==1){
        printf("uN^2 fm^%i\n",(2*L) - 2);
      }
    }
      
  }

  lt=lt*1.0E-12; //convert lifetime to s
  Et=Et/1000.0; //convert energy to MeV

  if(verbose){
    printf("\nB(%s) CALCULATION\n-----------------\n",mstr);
  }

  double fac = 8.0*PI*(L+1)/(L*HBAR_MEVS*pow(dblfac((2.0*L)+1.0),2.0)*LN2);
  fac = fac * (pow((Et/HBARC_MEVFM),2.0*L + 1.0));
  if(calcMode==0){
    /* lifetime to reduced transition probability */
    if(EM==0){
      /* electric */
      /* e^2 fm^(2L) */
      b=1/(fac*lt);
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
      b=1/(fac*lt);
      if(bup){
        b = b*(2.0*ji + 1.0)/(2.0*jf + 1.0);
      }
      if(L>1)
        printf("%0.4E uN^2 fm^%i\n",b,(2*L) - 2);
      else
        printf("%0.4E uN^2\n",b);
    }
  }else if(calcMode == 1){
    /* reduced transition probability to lifetime */
    if(EM==1){
      fac = fac * UN_MEVFM3/ESQ_MEVFM;
    }
    lt = 1/(fac*b);
    lt = lt / 1.0E-12; //convert lifetime to ps
    printf("%0.4E ps\n",lt);
  }

  return 0;
}