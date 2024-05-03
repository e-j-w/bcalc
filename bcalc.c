#include "bcalc.h"

void printHelp(){
  printf("\nReduced transition probability calculator\n");
  printf("example usage: bcalc -e VALUE -m VALUE -lt VALUE\n\n");
  printf("  Both of the following are needed:\n");
  printf("    -e         --  transition energy in keV\n");
  printf("    -m         --  multipole (eg. E1, M1, E2, etc.)\n");
  printf("\n");
  printf("  One of the following is needed:\n");
  printf("    -lt        --  Mean transition lifetime (in ps, use\n");
  printf("                   -ltns / -ltus / -lts / -lth for\n");
  printf("                   nanoseconds / microseconds / seconds / hours,\n");
  printf("                   respectively).\n");
  printf("    -hl        --  Transition half-life (in ps, use\n");
  printf("                   -hlns / -hlus / -hls / -hlh for\n");
  printf("                   nanoseconds / microseconds / seconds / hours,\n");
  printf("                   respectively).\n");
  printf("    -b         --  Reduced transition probability (for the\n"); 
  printf("                   L multipole) in units of e^2 fm^(2L) for\n");
  printf("                   electric multipoles or uN^2 fm^(2L-2) for\n");
  printf("                   magnetic multipoles.\n");
  printf("\n");
  printf(" --- Press any key for more ---");
  getc(stdin);
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -br        --  branching fraction of this transition (maximum 1,\n"); 
  printf("                   default 1)\n");
  printf("    -d         --  mixing ratio with the L+1 multipole\n");
  printf("    -icc       --  internal conversion coefficient for this transition\n");
  printf("                   (default 0)\n");
  printf("    -ji        --  inital spin (integer or half-integer)\n");
  printf("    -jf        --  final spin (integer or half-integer)\n");
  printf("    -A         --  mass number of the nucleus\n");
  printf("    -Z         --  proton number of the nucleus\n");
  printf("\n");
  printf(" --- Press any key for more ---");
  getc(stdin);
  printf("\n");
  printf("  Flags:\n");
  printf("    --barn     --  Use/calculate transition probability with spatial\n");
  printf("                   dimension in barns rather than fm (eg. e^2 b^2\n");
  printf("                   rather than e^2 fm^4).\n");
  printf("    --wu       --  Use/calculate transition probability in\n");
  printf("                   Weisskopf units (W.u.) rather than the default\n"); 
  printf("                   units specified above.  If used, requires the -A\n");
  printf("                   parameter.\n");
  printf("    --up       --  Use/calculate transition probability from\n");
  printf("                   final to initial state instead of vice versa.\n");
  printf("                   If used, requires the -ji and -jf parameters.\n");
  printf("    --brrel    --  Specifies that the branching fraction provided\n");
  printf("                   with the -br option is actually an intensity\n");
  printf("                   relative to another transition.\n");
  printf("    --beta2    --  Calculate the quadrupole deformation parameter,\n");
  printf("                   assuming a 2->0 (g.s.) E2 transition.  Requires\n");
  printf("                   '-m E2 -ji 2 -jf 0', and the -A and -Z parameters.\n");
  printf("                   Assumes mean charge radius R = r_0*A^(1/3), with\n");
  printf("                   r_0 = 1.2 fm.\n");
  printf("    --quiet    --  Only show the result of the calculation.\n");
}

double dblfac(unsigned int n){ 
  double val = 1.0;
  int i;
  for (i=(int)n; i>=0; i=i-2){
    if (i==0 || i==1)
      return val;
    else
      val *= i;
  }
  return val;
}

/* calculates single particle lifetimes */
double ltsp(const int EM, const int L, const int nucA, const double Et_keV){

  double hl_sp = 0.;
  double hbar = 6.58212E-19; /* keV s */
  double hbarc = 197.327E-10; /* keV cm */
  double esq = 1.440E-10; /* kev cm */
  double muNsq = 1.5922E-38; /* keV cm^3 */
  double RFac = 1.2E-13; /* cm */
  if(EM == 0){
    /* electric */
    hl_sp = log(2.0)*L*pow(dblfac((unsigned int)(2*L + 1)),2.0)*hbar*pow(((3.0+L)/3.0),2.0)*pow(hbarc,2*L + 1)/(2.0*(L + 1.0)*esq*pow(RFac,2.0*L));
    hl_sp /= (pow(Et_keV,2.0*L + 1)*pow(nucA,2.0*L/3.0));
  }else{
    /* magnetic */
    hl_sp = log(2.0)*L*pow(dblfac((unsigned int)(2*L + 1)),2.0)*hbar*pow(((3.0+L)/3.0),2.0)*pow(hbarc,2*L + 1)/(80.0*(L + 1.0)*muNsq*pow(RFac,2.0*L - 2.0));
    hl_sp /= (pow(Et_keV,2.0*L + 1)*pow(nucA,(2.0*L - 2.0)/3.0));
  }

  //convert from half-life to lifetime (in s)
  return hl_sp/LN2;
}

/* calculates the value of the quadrupole deformation parameter, assuming an input reduced transtion probability and a 2->0 transition */
void calcBeta2(const double Et, const double b_in, const int nucA, const int nucZ, const int barn, const int verbose){
  
  double b=0.;

  if(barn == 2){
    /* input using Weisskopf units, first calculate lifetime */
    double lt = ltsp(0,2,nucA,Et*1000.)/b_in;
    /* calculate b from lifetime */
    double fac = 8.0*PI*(2+1)/(2*HBAR_MEVS*pow(dblfac((2.0*2)+1.0),2.0)*LN2);
    fac = fac * (pow((Et/HBARC_MEVFM),2.0*2 + 1.0));
    b=1/(fac*lt); /* lifetime to reduced transition probability (E2: e^2 fm^4) */
  }else{
    b=b_in;
  }

  double beta = sqrt(5*b)*4.0*PI/(2*nucZ*ESQ_MEVFM*1.20*1.20*pow(1.0*nucA,2.0/3.0));

  if(verbose){
    printf("\nbeta_2 CALCULATION\n-----------------\n");
    printf("%0.4E\n",beta);
  }else{
    printf("beta_2 = %0.4E\n",beta);
  }
  
}

/* calculates the value of the quadrupole deformation parameter, assuming an input lifetime and a 2->0 transition */
void calcBeta2Lt(const double Et, const double lt, const int nucA, const int nucZ, const int verbose){

  double fac = 8.0*PI*(2+1)/(2*HBAR_MEVS*pow(dblfac((unsigned int)((2.0*2)+1.0)),2.0)*LN2);
  fac = fac * (pow((Et/HBARC_MEVFM),2.0*2 + 1.0));
  
  double b=1/(fac*lt); /* lifetime to reduced transition probability (E2: e^2 fm^4) */
  calcBeta2(Et,b,nucA,nucZ,0,verbose);
  
}

/* calculates the value of the reduced transtion probability */
void calcB(const int bup, const int EM, const int L, const double Et, const double lt, const double ji, const double jf, const int verbose, const int barn, const char *mstr, const int nucA){
  if(verbose){
    printf("\nB(%s) CALCULATION\n-----------------\n",mstr);
  }

  if(barn == 2){
    /* Weisskopf unit calculation */
    double lt_sp = ltsp(EM,L,nucA,Et*1000.);
    printf("%0.4E W.u.\n",lt_sp/lt);
    return;
  }

  double fac = 8.0*PI*(L+1)/(L*HBAR_MEVS*pow(dblfac((unsigned int)((2.0*L)+1.0)),2.0)*LN2);
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
      switch(barn){
        case 1:
          printf("%0.4E e^2 b^%i\n",b/pow(BARN_FM,L),L);
          break;
        case 0:
        default:
          printf("%0.4E e^2 fm^%i\n",b,2*L);
          break;
      }  
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
    if(L>1){
      switch(barn){
        case 1:
          printf("%0.4E uN^2 b^%i\n",b/pow(BARN_FM,L-1),L-1);
          break;
        case 0:
        default:
          printf("%0.4E uN^2 fm^%i\n",b,(2*L) - 2);
          break;
      }
    }else{
      printf("%0.4E uN^2\n",b);
    }
      
  }
}

void calcLt(const int bup, const int EM, const int L, const double Et, double b, const double ji, const double jf, const int verbose, const int barn, const char *mstr, const int nucA, const double branching){
  if(verbose){
    printf("\nLIFETIME CALCULATION\n--------------------\n");
  }

  if(bup){
    b = b*(2.0*jf + 1.0)/(2.0*ji + 1.0);
  }

  if(barn == 2){
    /* Weisskopf unit calculation */
    double lt_sp = ltsp(EM,L,nucA,Et*1000.);
    lt_sp = lt_sp / 1.0E-12; //convert from s to ps
    printf("%0.4E ps\n",lt_sp/b);
    return;
  }

  //convert barn units to fm units
  if(barn){
    if(EM==0){
      b = b*pow(BARN_FM,L);
    }else if(EM==1){
      b = b*pow(BARN_FM,L-1);
    }
  }

  double fac = 8.0*PI*(L+1)/(L*HBAR_MEVS*pow(dblfac((unsigned int)((2.0*L)+1.0)),2.0)*LN2);
  fac = fac * (pow((Et/HBARC_MEVFM),2.0*L + 1.0));
  /* reduced transition probability to lifetime */
  if(EM==1){
    fac = fac * UN_MEVFM3/ESQ_MEVFM;
  }
  double lt = 1/(fac*b);
  lt = lt / 1.0E-12; //convert lifetime to ps
  lt = 1.0/((1.0/lt)*branching); //partial lifetime
  printf("%0.4E ps",lt);
  if(branching != 1.){
    printf(" (partial lifetime)");
  }
  printf("\n");
}

int main(int argc, char *argv[]) {

  if (argc == 1) {
    printHelp();
    exit(-1);
  }

  int i; /*counters*/

  /*initialize parameter values*/
  char mstr[4], mstr1[12];
  double Et = -1.; /* transition energy */
  int L = -1; /* multipolarity */
  int EM = -1; /* 0=electric, 1=magnetic */
  double lt = 0.; /* lifetime */
  double lt1 = 0.; /* lifetime for L+1 multipole (mixed transitions) */
  double b = 0.; /* reduced transition probability */
  int calcMode = -1; /* 0=calulate B, 1=calculate lifetime */
  int verbose = 1; /* 0=none, 1=verbose */
  int barn = 0; /* 0=fm units, 1=barn units, 2=Weisskopf units */
  int bup = 0; /* 0=down, 1=up */
  double ji = -1.; /* Initial spin */
  double jf = -1.; /* Final spin*/
  double delta = 0; /* mixing ratio */
  double branching = 1.; /* branching fraction */
  double icc = 0.; /* internal conversion coefficient */
  int useDelta = 0; /* 0=no mixing, 1=mixing */
  int calcB2 = 0; /* 0=no beta_2 calc, 1=calc beta_2 */
  int brrel = 0; /* 0=use branching fraction, 1=use relative intensity */
  int nucA = -1; /* mass number of the nucleus of interest */
  int nucZ = -1; /* proton number of the nucleus of interest */

  /*read parameters*/
  for(i=0;i<argc;i++){
    if(strcmp(argv[i],"--quiet")==0){
      verbose = 0;
    }else if(strcmp(argv[i],"--up")==0){
      bup = 1;
    }else if(strcmp(argv[i],"--beta2")==0){
      calcB2 = 1;
    }else if(strcmp(argv[i],"--barn")==0){
      barn += 1;
    }else if(strcmp(argv[i],"--wu")==0){
      barn += 2;
    }else if(strcmp(argv[i],"--brrel")==0){
      brrel = 1;
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
      mstr[2] = '\0';
      strncpy(mstr,argv[i+1],sizeof(mstr));
      if(mstr[0] == 'E'){
        EM=0;
      }else if(mstr[0] == 'M'){
        EM=1;
      }else{
        printf("ERROR: invalid multipole value.\n");
        exit(-1);
      }
      if((isdigit(mstr[2])!=0)&&(isdigit(mstr[1])!=0)){
        L=mstr[2] - '0';
        L+=(mstr[1] - '0')*10;
      }else if(isdigit(mstr[1])!=0){
        L=mstr[1] - '0';
      }else{
        printf("ERROR: invalid multipole value.\n");
        exit(-1);
      }
      
    }else if((strcmp(argv[i],"-Lt")==0)||(strcmp(argv[i],"-lt")==0)||(strcmp(argv[i],"-Ltps")==0)||(strcmp(argv[i],"-ltps")==0)){
      lt=atof(argv[i+1]);
      calcMode = 0;
    }else if((strcmp(argv[i],"-Ltns")==0)||(strcmp(argv[i],"-ltns")==0)){
      lt=atof(argv[i+1])*1000.0;
      calcMode = 0;
    }else if((strcmp(argv[i],"-Ltus")==0)||(strcmp(argv[i],"-ltus")==0)){
      lt=atof(argv[i+1])*1000000.0;
      calcMode = 0;
    }else if((strcmp(argv[i],"-Lts")==0)||(strcmp(argv[i],"-lts")==0)){
      lt=atof(argv[i+1])*1000000000000.0;
      calcMode = 0;
    }else if((strcmp(argv[i],"-Lth")==0)||(strcmp(argv[i],"-lth")==0)){
      lt=atof(argv[i+1])*3600*1000000000000.0;
      calcMode = 0;
    }else if((strcmp(argv[i],"-Hl")==0)||(strcmp(argv[i],"-hl")==0)||(strcmp(argv[i],"-Hlps")==0)||(strcmp(argv[i],"-hlps")==0)){
      lt=atof(argv[i+1])/LN2;
      calcMode = 0;
    }else if((strcmp(argv[i],"-Hlns")==0)||(strcmp(argv[i],"-hlns")==0)){
      lt=atof(argv[i+1])*1000.0/LN2;
      calcMode = 0;
    }else if((strcmp(argv[i],"-Hlus")==0)||(strcmp(argv[i],"-hlus")==0)){
      lt=atof(argv[i+1])*1000000.0/LN2;
      calcMode = 0;
    }else if((strcmp(argv[i],"-Hls")==0)||(strcmp(argv[i],"-hls")==0)){
      lt=atof(argv[i+1])*1000000000000.0/LN2;
      calcMode = 0;
    }else if((strcmp(argv[i],"-Hlh")==0)||(strcmp(argv[i],"-hlh")==0)){
      lt=atof(argv[i+1])*3600*1000000000000.0/LN2;
      calcMode = 0;
    }else if((strcmp(argv[i],"-B")==0)||(strcmp(argv[i],"-b")==0)){
      b=atof(argv[i+1]);
      calcMode = 1;
    }else if(strcmp(argv[i],"-d")==0){
      delta=atof(argv[i+1]);
      useDelta = 1;
    }else if(strcmp(argv[i],"-br")==0){
      branching=atof(argv[i+1]);
    }else if(strcmp(argv[i],"-icc")==0){
      icc=atof(argv[i+1]);
      if(icc < 0.){
        printf("ERROR: Internal conversion coefficient must be a positive number.\n");
        exit(-1);
      }
    }else if(strcmp(argv[i],"-A")==0){
      nucA=atoi(argv[i+1]);
    }else if(strcmp(argv[i],"-Z")==0){
      nucZ=atoi(argv[i+1]);
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
  if(L > 12){
    printf("ERROR: Maximum calculable L-value is 12.\n");
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
  if((branching <= 0.)||(branching > 1.)){
    printf("ERROR: invalid branching fraction (must be a positive number less than 1).\n");
    exit(-1);
  }
  if(barn>2){
    printf("ERROR: only one of --barn and --wu can be used at once.\n");
    exit(-1);
  }
  if(barn==2){
    if(nucA == -1){
      printf("ERROR: when using --wu, must also specify the -A parameter.\n");
      exit(-1);
    }else if(nucA <= 0){
      printf("ERROR: The mass number A cannot be less than 1.\n");
      exit(-1);
    }
  }
  if(calcB2){
    if((EM!=0)||(L!=2)||(ji!=2)||(jf!=0)){
      printf("ERROR: Can only calculate beta_2 for 2->0 (g.s.) E2 transitions.  The parameter values '-m E2 -ji 2 -jf 0' are required.\n");
      exit(-1);
    }
    if((nucZ == -1)||(nucA == -1)){
      printf("ERROR: when calculating beta_2, must specify both the -A and -Z parameters.\n");
      exit(-1);
    }else if((nucZ <= 0)||(nucA <= 0)){
      printf("ERROR: The proton number Z or mass number A cannot be less than 1.\n");
      exit(-1);
    }
  }
  if(nucA < nucZ){
    printf("ERROR: The mass number A cannot be less than the proton number Z.\n");
    exit(-1);
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
    else if(L==5)
      printf("triacontadipole");
    else if(L==6)
      printf("hexacontatetrapole");
    else if(L==7)
      printf("hecatonicosioctopole");
    else if(L==8)
      printf("diacosiapentecontahexadecapole");
    else
      printf("L = %i",L);
    if(useDelta&&(calcMode==0)){
      if(delta > 0.01){
        printf(" (L+1 mixing, delta = %0.3f)", delta);
      }else if(delta > 0.00001){
        printf(" (L+1 mixing, delta = %0.6f)", delta);
      }else{
        printf(" (L+1 mixing, delta = %0.9f)", delta);
      }
    }
    printf("\n");
    if(icc > 0.){
      printf("Internal conversion coefficient: %0.3f\n",icc);
    }
    if(calcMode == 0){
      if(lt < 1E3){
        printf("Mean lifetime: %0.3f ps\n",lt);
      }else if(lt < 1E12){
        printf("Mean lifetime: %0.3f ns\n",lt/((double)1E3));
      }else if(lt < 1E16){
        printf("Mean lifetime: %0.3f s\n",lt/((double)1E12));
      }else
        printf("Mean lifetime: %0.3f hr\n",lt/(3600.0*(double)1E12));
    }else if(calcMode == 1){
      printf("B(%s): %f ",mstr, b);
      if(EM==0){
        if(barn == 0){
          printf("e^2 fm^%i\n",2*L);
        }else if(barn == 1){
          printf("e^2 b^%i\n",L);
        }else{
          printf("W.u.\n");
        }
      }else if(EM==1){
        if(L>1){
          if(barn == 0){
            printf("uN^2 fm^%i\n",(2*L) - 2);
          }else if(barn == 1){
            printf("uN^2 b^%i\n",L-1);
          }else{
            printf("W.u.\n");
          }
        }else if(L==1){
          if(barn < 2){
            printf("uN^2\n");
          }else{
            printf("W.u.\n");
          }
        }
      }
    }
    if(nucA>0)
      printf("A = %i\n",nucA);
    if(nucZ>0)
      printf("Z = %i\n",nucZ);
    if(brrel == 0){
      if(branching > 0.1){
        printf("Branching fraction: %.2f\n",branching);
      }else if(branching > 0.001){
        printf("Branching fraction: %.4f\n",branching);
      }else if(branching > 0.0000001){
        printf("Branching fraction: %.8f\n",branching);
      }else if(branching > 0.00000000001){
        printf("Branching fraction: %.12f\n",branching);
      }else{
        printf("Branching fraction: %.16f\n",branching);
      }
    }else if(brrel == 1){
      if(branching > 0.1){
        printf("Relative intensity: %.2f\n",branching);
      }else if(branching > 0.001){
        printf("Relative intensity: %.4f\n",branching);
      }else if(branching > 0.0000001){
        printf("Relative intensity: %.8f\n",branching);
      }else if(branching > 0.00000000001){
        printf("Relative intensity: %.12f\n",branching);
      }else{
        printf("Relative intensity: %.16f\n",branching);
      }
    }
  }

  /* branching fraction calculation */
  if(brrel == 1)
    branching = branching/(branching + 1.0);
  if(calcMode == 0){
    lt = 1.0/((1.0/lt)*branching); //partial lifetime
    lt = lt*(1.0 + icc);
  }

  /* mixing ratio calculation */
  if(useDelta){
    lt1 = lt * (1.0 + delta*delta) / (delta*delta);
    lt = lt * (1.0 + delta*delta);
    if(mstr[0] == 'E')
      snprintf(mstr1,12,"M%i",L+1);
    else
      snprintf(mstr1,12,"E%i",L+1);
  }

  /* report partial lifetimes */
  if((verbose)&&(calcMode==0)){
    if((branching != 1.)&&(!useDelta)){
      if(lt < 1E3){
        printf("Partial lifetime: %0.3f ps\n",lt);
      }else if(lt < 1E12){
        printf("Partial lifetime: %0.3f ns\n",lt/((double)1E3));
      }else if(lt < 1E16){
        printf("Partial lifetime: %0.3f s\n",lt/((double)1E12));
      }else{
        printf("Partial lifetime: %0.3f hr\n",lt/(3600.0*(double)1E12));
      }
    }else if(useDelta){
      if(lt < 1E3){
        printf("Partial lifetime (%s): %0.3f ps\n",mstr,lt);
      }else if(lt < 1E12){
        printf("Partial lifetime (%s): %0.3f ns\n",mstr,lt/((double)1E3));
      }else if(lt < 1E16){
        printf("Partial lifetime (%s): %0.3f s\n",mstr,lt/((double)1E12));
      }else{
        printf("Partial lifetime (%s): %0.3f hr\n",mstr,lt/(3600.0*(double)1E12));
      }
      if(lt1 < 1E3){
        printf("Partial lifetime (%s): %0.3f ps\n",mstr1,lt1);
      }else if(lt1 < 1E12){
        printf("Partial lifetime (%s): %0.3f ns\n",mstr1,lt1/((double)1E3));
      }else if(lt1 < 1E16){
        printf("Partial lifetime (%s): %0.3f s\n",mstr1,lt1/((double)1E12));
      }else{
        printf("Partial lifetime (%s): %0.3f hr\n",mstr1,lt1/(3600.0*(double)1E12));
      }
    }
  }


  lt=lt*1.0E-12; //convert lifetime to s
  lt1=lt1*1.0E-12; //convert lifetime to s
  Et=Et/1000.0; //convert energy to MeV


  if(calcMode == 0){
    calcB(bup,EM,L,Et,lt,ji,jf,verbose,barn,mstr,nucA);
    if(useDelta){
      calcB(bup,!EM,L+1,Et,lt1,ji,jf,verbose,barn,mstr1,nucA);
    }
    if(calcB2){
      calcBeta2Lt(Et,lt,nucA,nucZ,verbose);
    }
  }else if(calcMode == 1){
    calcLt(bup,EM,L,Et,b,ji,jf,verbose,barn,mstr,nucA,branching);
    if(calcB2){
      calcBeta2(Et,b,nucA,nucZ,barn,verbose);
    }
  }
  

  return 0;
}