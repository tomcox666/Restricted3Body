/******************************************************************/
/*                                                                */
/* RESTRICTED 3-BODY PROBLEM INTEGRATOR                           */
/*                                                                */
/* RICHARD ALEXANDER                                              */
/* ORIGINALLY WRITTEN DEC 2017; THIS VERSION DEC 2019             */
/*                                                                */
/* INTEGRATES TEST PARTICLE ORBITING AROUND FIXED BINARY          */
/* LEAPFROG INTEGRATOR; KICK-DRIFT-KICK FORM (SYNCRHONISED)       */
/*                                                                */
/* FIXED GLOBAL TIMESTEPPING                                      */
/* OUTPUT DUMPS AT FIXED INTERVALS                                */
/*                                                                */
/* VARIABLES:                                                     */
/* X[i][] IS THE POSITION VECTOR OF THE ith STAR; M[i] IS MASS    */
/* x[],v[],a[] ARE POS, VEL, ACC OF PLANET (TEST PARTICLE)        */
/*                                                                */
/* UNITS ARE SCALE-FREE: M[0]+M[1] = 1, a_b=1, G=1                */
/* => BINARY ORBITAL PERIOD = 2pi (CODE UNITS)                    */
/*                                                                */
/* FORCE COMPUTATION SUB-ROUTINE INCLUDES POTENTIAL ENERGY        */
/* OUTPUTS ORBITAL ELEMENTS: E & J NEEDED TO COMPUTE a & e        */
/* (e AT PRESENT ASSUMES CO-PLANAR ORBITS)                        */
/*                                                                */
/* BINARY HAS FIXED ECCENTRICITY e_b                              */
/* PLANET HAS RANDOM INITIAL PHASE ANGLE                          */
/* COMPUTES a & e FOR PLANET (ASSUMING KEP ORBIT & POTENTIAL)     */
/*                                                                */
/* TAKES BINARY PERIOD IN DAYS AS AN INPUT, SO ON-SCREEN OUTPUT   */
/* IS IN PHYSICAL UNITS. BUT THIS IS ONLY FOR SCREEN OUTPUT.      */
/* INTERNAL CALCULATIONS AND ALL OUTPUT DUMPS ARE IN CODE UNITS.  */
/*                                                                */
/*                                                                */
/*                                                                */
/*                                                                */
/* COMPILE WITH gcc, SYNTAX:                                      */
/*           gcc restricted_3body.c -o restricted_3body.e -lm     */
/* THEN RUN:                                                      */
/*     ./restricted_3body.e -p P -q Q -e E -r R -t T -d D -l L    */
/*                                                                */
/* REQUIRED INPUT PARAMETERS ARE:                                 */
/*       Q = binary mass ratio                                    */
/*       E = binary eccentricity                                  */
/*       R = planet semi-major axis (units of binary separation)  */
/* OPTIONAL INPUT PARAMETERS:                                     */
/*       P = binary period (days; default = 10)                   */
/*       T = no of timesteps per binary orbit (default = 100)     */
/*       D = timesteps between output dumps (default = 10^4)      */
/*       L = time limit in binary orbits (default = 10^6)         */
/*                                                                */
/*                                                                */
/* EXAMPLE: ./restricted_3body.e -q 1.0 -e 0.0 -r 3.0             */
/* THIS WILL SET UP AN EQUAL-MASS (q=1) CIRCULAR (e=0) BINARY,    */
/* WITH THE PLANET AT 3 TIMES THE BINARY SEPARATION, AND THE      */
/* OTHER (OPTIONAL) PARAMETERS SET TO THE DEFAULTS.               */
/*                                                                */
/*                                                                */
/*                                                                */
/* OUTPUT FILES ARE IN THE data/ SUB-DIRECTORY                    */
/* ************************************************************** */
/* **** IF THIS DIRECTORY DOES NOT EXIST THE CODE WILL CRASH **** */
/* **** ENTER COMMAND mkdir data TO CREATE IT (IF NEEDED) ******* */    
/* ************************************************************** */
/*                                                                */
/* OUTPUT FILES:                                                  */
/*  data/planet_*.dat - pos/vel/acc/orbit constants of planet     */
/*  data/binary_*.dat - masses & positions of both stars          */
/*  data/*.units      - log file with optional physical units     */
/*                                                                */
/* FILE-NAMES ARE AUTOMATICALLY GENERATED FROM "HARD-CODED" ROOT  */
/* STRING AND THE INPUT PARAMETERS (READ FROM COMMAND LINE).      */
/*                                                                */
/* TIME-STEPPING - DEFAULT IS dt = P_b/100, BUT THIS SHOULD BE    */
/*                 TESTED ON A CASE-BY-CASE BASIS.                */
/*                                                                */
/* CAN CALCULATE UP TO ~1E7 TIMESTEPS PER CPU-SECOND, BUT REAL    */
/* RUN-TIMES DEPEND ON OUTPUT FREQUENCY.                          */
/* HIGH OUTPUT FREQUENCY WILL SLOW CODE DOWN BY A LARGE FACTOR.   */
/*                                                                */
/******************************************************************/

/* STANDARD HEADER FILES */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/* FORCE COMPUTATION SUBROUTINE */
void compute_forces(int N, double X[2][3], double M[2], double x[3], double a[3], double*U)
{
  int i,j;
  double r,r2,r3;
  double UU;

  /* SET ACCELERATION VECTOR TO ZERO */
  for (j=0; j<3; j++)
    a[j] = 0.0;
  /* SET POTENTIAL ENERGY TO ZERO */
  UU = 0.0;

  /* LOOP OVER STAR PARTICLES */
  for (i=0; i<N; i++)
    {
      /* COMPUTE r^3 FOR STAR i */
      r2 = 0.0;      
      for (j=0; j<3; j++)
	r2 = r2 + (x[j]-X[i][j])*(x[j]-X[i][j]);
      r = sqrt(r2);
      r3 = r2 * r;

      /* INCREMENT POTENTIAL ENERGY */
      UU = UU - M[i]/r;

      /* INCREMENT ACCELERATIONS FOR STAR i */
      for (j=0; j<3; j++)
	a[j] = a[j] - (M[i]/r3)*(x[j]-X[i][j]);

    }
  /* POINTER FOR POTENTIAL ENERGY */
  *U = UU;

  return;

}


/* RANDOM NUMBER GENERATOR */
double rand_num(double min, double max)
{
    double random_num;
       
    random_num = ( (double)rand() / ((double)(RAND_MAX)) );   
   
    random_num = (random_num * (max-min)) + min;   // output # from [r_min, r_max]
   
    return (random_num);
}


/* MAIN PROGRAM */
int main(int argc, char** argv)
{
  /* DECLARE ALL VARIABLES */
  int i,j,N,dt_dump,steps_to_dump,unstable;
  double Time,dt,time_limit;
  /* STARS */
  double X[2][3];
  double Omega_b;
  double M[2];
  /* PLANET */
  double x[3];
  double v[3];
  double a[3];
  double m;
  /* OTHER VARIABLES */
  double r2,r3;
  double Mtot,a_b,q,a0,a1,P,e_b;
  double U,K,E,E_tot;
  double a_p,e_p,r;
  double delta,v_circ,r_p,J_z,e2;
  double time_to_dump;


  /* DEFINE CONSTANTS IN cgs UNITS */
  static const double G_cgs = 6.672E-8;
  static const double Msun = 1.989E33;
  static const double Mjup = 1.898E30;
  static const double AU = 1.496E13;
  static const double Rsun = 6.955E10;
  static const double Rjup = 6.9911E9;
  static const double year = 3.15569E7;
  static const double day = 86400.0;
  static const double sigma_p = 5.9E-54;
  static const double sigma_s = 6.4E-59;
  /* VARIABLES IN PHYSICAL UNITS */
  double P_binary,P_cgs,a_b_cgs,Omega_cgs;


  /* COMMAND LINE INPUT VARIABLES */
  float massrat,eccen,P_b,period,r_planet,n_steps,n_orbits;
  int out_steps;


  /* OUTPUT FILES, ETC. */
  char datafile[100],unitsfile[100],root[100];
  FILE *planet_data,*binary_data,*units;

  /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
  /* PREAMBLE: READ PARAMETERS & OPEN OUTPUT FILES */
  /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/


  /* SET DEFAULTS FOR OPTIONAL INPUT PARAMETERS */
  period = 10.0;        /* BINARY PERIOD IN DAYS */
  n_steps = 50.0;     /* NO OF TIMESTEPS PER BINARY ORBIT */
  out_steps = 500.0;   /* NO OF TIMESTEPS BETWEEN OUTPUT DUMPS */
  n_orbits = 1.0E6;    /* TIME LIMIT (STOPPING POINT), IN BINARY ORBITS */
  
  /* READ INPUT PARAMETERS FROM COMMAND LINE */
  for (i=1; i<argc; i++)
    {

      if (argv[i][0] == '-')
	{
	  switch (argv[i][1])
	    {
	      /* REQUIRED INPUTS */
	    case 'q':   massrat = atof(argv[++i]);   /* BINARY MASS RATIO */
	      break;
	    case 'e':   eccen = atof(argv[++i]);     /* BINARY ECCENTRICITY */
	      break;
	    case 'r':   r_planet = atof(argv[++i]);  /* PLANET SEMI-MAJ AXIS */
	      break;                                 /* (UNITS OF a_bin)     */
	      /* OPTIONAL INPUTS */
	    case 'p':   period = atof(argv[++i]);    /* BINARY PERIOD IN DAYS */
	      break;
	    case 't':   n_steps = atof(argv[++i]);   /* TIMESTEPS PER ORBIT */
	      break;
	    case 'd':   out_steps = atoi(argv[++i]); /* TIMESTEPS PER OUTPUT */
	      break;
	    case 'l':   n_orbits = atof(argv[++i]);  /* TIME LIMIT */
	      break;
	      
	    }
	}
    }

  /* CREATE OUTPUT FILENAME */
  sprintf(root,"res3body_e=%.5f",eccen);
  printf("%s\n",root);

  /* OPEN OUTPUT FILES */
  sprintf(datafile,"data/planet_%s.txt",root);
  planet_data = fopen(datafile,"w");
  sprintf(datafile,"data/binary_%s.txt",root);
  binary_data = fopen(datafile,"w");

  /* INITIALISE RANDOM NUMBER GENERATOR WITH TIME OF DAY */
  srand(time(NULL));
  delta = rand();


  /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

  /* STEP 0 - SET INITIAL CONDITIONS */

  /* STARS */
  /* NUMBER OF STARS */
  N = 2;
  /* TOTAL MASS */
  Mtot = 1.0;
  /* BINARY SEPARATION */
  a_b = 1.0;

  /* SET PERIOD, MASS RATIO AND ECCENTRICITY (READ FROM COMMAND LINE) */
  q = massrat;
  e_b = eccen;
  P_b = period;

  /* DEFINE BINARY ORBITAL FREQUENCY */
  /* STRICTLY THIS IS UNNECESSARY, AS Omega_b==1, BUT GOOD PRACTICE */
  Omega_b = sqrt(Mtot/(a_b*a_b*a_b));

    
  /* SET N-BODY ICs*/
  /* BINARY */
  /* COMPUTE STELLAR MASSES FROM MASS RATIO */
  M[0] = Mtot / (1.0+q) ;
  M[1] = Mtot - M[0];
  /* FIND DISTANCES FROM CENTRE OF MASS */
  a1 = a_b / (1.0+q);
  a0 = a_b - a1;

  /* SET TIME TO ZERO */
  Time = 0.0;
  /* SET (INTIAL POSITIONS) OF THE STARS */
  X[0][0] = -1.0 * a0 * ((1.0-e_b*e_b)/(1.0+e_b*cos(Omega_b*Time))) * cos(Omega_b*Time);
  X[0][1] = -1.0 * a0 * ((1.0-e_b*e_b)/(1.0+e_b*cos(Omega_b*Time))) * sin(Omega_b*Time);
  X[0][2] = 0.0;
  
  X[1][0] = a1 * ((1.0-e_b*e_b)/(1.0+e_b*cos(Omega_b*Time))) * cos(Omega_b*Time);
  X[1][1] = a1 * ((1.0-e_b*e_b)/(1.0+e_b*cos(Omega_b*Time))) * sin(Omega_b*Time);
  X[1][2] = 0.0;
  
  
  /* PLANET - KEPLERIAN CIRCULAR VELOCITY AT r_p */
  /* ASSIGN RANDOM INITIAL ORBITAL PHASE */
  delta = rand_num(0,1) * 2.0 * M_PI;
  printf("phi_0 = %g\n",delta);

  /* SET PLANET SEMI-MAJOR AXIS (READ FROM COMMAND LINE) */
  r_p = r_planet;
  /* CALCULATE CIRCULAR (KEPLERIAN) VELOCITY */
  v_circ = sqrt(((M[0]+M[1]))/r_p);
  printf("r_p = %g, v_c = %g\n",r_p,v_circ);

  /* SET PLANET ICS */
  x[0] = r_p * cos(delta);
  x[1] = r_p * sin(delta);
  x[2] = 0.0;
  v[0] = -1.0 * v_circ * sin(delta);
  v[1] = v_circ * cos(delta);
  v[2] = 0.0;
  /* COMPUTE DISTANCE FROM CENTRE OF MASS */
  r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);


  /* SET TIMESTEP (AS FRACTION OF BINARY PERIOD) */
  P = 2.0*M_PI / Omega_b;
  dt = P / n_steps;
  /* SET OUTPUT FREQUENCY (EVERY dt_dump TIMESTEPS) */
  dt_dump = out_steps;

  /* INITIALISE INTEGRATOR */
  Time = 0.0;
  unstable = 0;
  /* INITIAL FORCE COMPUTATION */
  compute_forces(N, X, M, x, a, &U);

  /* COMPUTE ENERGY OF PLANET (IN JACOBI COORDINATES) */
  K = 0.5*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  U = -1.0*(M[0]+M[1])/r;
  E = U + K;

   /* FIND SEMI-MAJOR AXIS */
  a_p = -0.5*(M[0]+M[1])/E;
  /* FIND ECCENTRICITY */
  /* ANG MOM ASSUMES 2D, ONLY Z-COMPONENT */
  J_z = (x[0]*v[1] - x[1]*v[0]);
  e2 = 1.0 + 2.0*E*J_z*J_z/(M[0]+M[1]);
  if (e2 < 0)
    e_p = 0.0;
  else
    e_p = sqrt(e2);

  printf("a(0) = %g; e(0) = %g\n",a_p,e_p);

  /* DEFINE PHYSICAL UNITS - SPECIFY BINARY PERIOD */
  P_binary = P_b;    /* days */
  P_cgs = P_binary * day;   /* s */
  Omega_cgs = 2.0*M_PI/P_cgs;         /* Hz */
  a_b_cgs = pow( (G_cgs * Msun * Mtot)/(Omega_cgs*Omega_cgs) , (1.0/3.0) );
  printf("Binary orbital period = %gdays\n",P_binary);
  printf("Binary semi-major axis = %gAU\n",a_b_cgs/AU);
  printf("dt = %g P_b\n",dt/P);
  printf("dt_output = %g P_b\n",dt/P*dt_dump);

  /* SET TIME LIMIT - CODE WILL STOP AFTER n_orbits BINARY ORBITS */
  time_limit = 2.0*M_PI*n_orbits;
  printf("Time limit = %g (CODE UNITS); %g yr\n",time_limit,time_limit/(2.0*M_PI)*P_cgs/(year) );

  /* INITIAL OUTPUT DUMP */
  /* PLANET */
  fprintf(planet_data,"%g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g\n",Time,x[0],x[1],x[2],v[0],v[1],v[2],a[0],a[1],a[2],E,a_p,e_p);
  /* BINARY */
  fprintf(binary_data,"%g  %g  %g  %g  %g  %g  %g  %g\n",M[0],X[0][0],X[0][1],X[0][2],M[1],X[1][0],X[1][1],X[1][2]);


  /* CREATE UNITS FILE */
  sprintf(unitsfile,"data/%s.units",root);  
  units = fopen(unitsfile,"w");
  fprintf(units,"%g\n",(P_cgs/year)/(2.0*M_PI));
  fprintf(units,"%g\n",a_b_cgs);
  fclose(units);

  /* RESET OUTPUT COUNTER */
  steps_to_dump = dt_dump;

  /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/


  /* START TIME INTEGRATION */
  while ((Time <= time_limit) && (unstable == 0))
    {



      /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

      /* STEP 1 - KICK  */

      /* KICK - v(t) -> v(t+1/2) */
      for (j=0; j<3; j++)
	v[j] = v[j] + 0.5*dt*a[j];
      

      /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

      /* STEP 2 - DRIFT  */
      
      /* DRIFT - r(t) -> r(t+1) */
      for (j=0; j<3; j++)
	x[j] = x[j] + v[j]*dt;

      /* COMPUTE DISTANCE FROM CENTRE OF MASS */
      r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

      /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
      
      /* STEP 3 - MOVE BINARY  */
 
      /* INCREMENT TIME - d -> d+dt */
      Time = Time + dt;
      steps_to_dump--;
      
       /* MOVE STAR PARTICLES (ON FIXED ORBIT) - X(t) -> X(t+dt) */
      X[0][0] = -1.0 * a0 * ((1.0-e_b*e_b)/(1.0+e_b*cos(Omega_b*Time))) * cos(Omega_b*Time);
      X[0][1] = -1.0 * a0 * ((1.0-e_b*e_b)/(1.0+e_b*cos(Omega_b*Time))) * sin(Omega_b*Time);
      X[0][2] = 0.0;
      
      X[1][0] = a1 * ((1.0-e_b*e_b)/(1.0+e_b*cos(Omega_b*Time))) * cos(Omega_b*Time);
      X[1][1] = a1 * ((1.0-e_b*e_b)/(1.0+e_b*cos(Omega_b*Time))) * sin(Omega_b*Time);
      X[1][2] = 0.0;      


      /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
      
      /* STEP 4 - COMPUTE FORCES  */
      
      /* FORCE COMPUTATION */
      compute_forces(N, X, M, x, a, &U);

      /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

      /* STEP 5 - KICK  */
     
      /* KICK - v(t+1/2) -> v(t+1) */
      for (j=0; j<3; j++)
	v[j] = v[j] + 0.5*dt*a[j];
      
      /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

      /* STEP 6 - BOOK-KEEPING  */
      
      /* COMPUTE KINETIC ENERGY OF PLANET */
      K = 0.5*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
      /* COMPUTE TOTAL ENERGY */
      E_tot = U+K;
      /* FOR CALCULATING e & a, ASSUME POINT MASS POTENTIAL (JACOBI COORDS) */
      U = -1.0*(M[0]+M[1])/r;
      /* COMPUTE JACOBI TOTAL ENERGY */
      E = U + K;
      
      /* FIND SEMI-MAJOR AXIS */
      a_p = -0.5*(M[0]+M[1])/E;
      /* FIND ECCENTRICITY */
      /* ANG MOM ASSUMES 2D, ONLY Z-COMPONENT */
      /* (CAN GENERALISE TO 3D FOR OUT-OF-PLANE ORBITS) */
      J_z = (x[0]*v[1] - x[1]*v[0]);
      e2 = 1.0 + 2.0*E*J_z*J_z/(M[0]+M[1]);
      if (e2 < 0)
	e_p = 0.0;
      else
	e_p = sqrt(e2);


      /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/


      /* TEST FOR AUTOMATIC STOPPING CRITERIA        */
      /* CODE STOPS IF 1) PLANET INSIDE BINARY ORBIT */
      /*               2) PLANET ORBIT IS UNBOUND    */
      if ((r < a_b) || ((r > 50.0*a_b) && (E_tot > 0.0)) || (a_p/a_b > 50.0) || (e_p > 0.95) )
	unstable = 1;

      /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/


      /* STEP 7 - DATA OUTPUT  */

      /* OUTPUT STEP */
      if ((steps_to_dump <= 0) || (unstable == 1))
	{
	  /* OUTPUT PARTICLE DATA TO FILE */
	  /* PLANET */
	  fprintf(planet_data,"%g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g\n",Time,x[0],x[1],x[2],v[0],v[1],v[2],a[0],a[1],a[2],E,a_p,e_p);
	  /* BINARY */
	  fprintf(binary_data,"%g  %g  %g  %g  %g  %g  %g  %g\n",M[0],X[0][0],X[0][1],X[0][2],M[1],X[1][0],X[1][1],X[1][2]);
	  fflush(planet_data);
	  fflush(binary_data);

	  /* RESET COUNTER */
	  steps_to_dump = dt_dump;

	  /* SCREEN OUTPUT */
	  printf("t = %g yr  a_p = %g  e_p = %g\n",Time/(2.0*M_PI)*P_cgs/(year),a_p,e_p);
	}


      /*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
      
      
    }
  
  /* CLOSE OUTPUT DATA FILE */
  fclose(planet_data);
  fclose(binary_data);


  /* ECHO END STATE TO SCREEN */
  if (unstable == 0)
    printf("Stopped at time limit: t = %g yr\n",Time/(2.0*M_PI)*P_cgs/(year));
  if (unstable == 1)
    printf("Unstable: stopped at t = %g yr\n",Time/(2.0*M_PI)*P_cgs/(year));
  
  /* ECHO TO UNITS FILE */
  units = fopen(unitsfile,"a");
  fprintf(units,"%g  %i\n",Time,unstable);
  fclose(units);


  exit(EXIT_SUCCESS);
}

