/********************************************
 *
 * This code is the merge of mcxyz.c and the Appendix from the thesis of Zhao
 *
 * myname = 'skinvessel_mergedcode_2layers';
 *
 * The program reads two files prepared by user:
 *  myname_H.mci    = header input file for mcxyz
 *  myname_T.bin    = tissue structure file
 *
 * The output will be written to the following files:
 *
 *  myname_DetS.bin = path legths of detected photons
 *  myname_DetW.bin = weight of detected photons
 *  myname_DetL.bin = likelihood ratio of detected photons
 *  myname_DetID.bin = ID of detected photons
 *  myname_DetZ.bin = maximal depth of detected photons
 *  myname_F.bin    = fluence rate output F[i] [W/cm^2 per W delivered]
 *  myname_Ryx.bin = reflectance/escaping flux R[i] [W/cm^2 per W delivered]
 *  myname_Rd.dat =
 *
 **********/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define Ntiss		19          /* Number of tissue types. */
#define STRLEN 		32          /* String length. */
#define ls          1.0E-7      /* Moving photon a little bit off the voxel face */
#define	PI          3.1415926
#define	LIGHTSPEED	2.997925E10 /* in vacuo speed of light [cm/s] */
#define ALIVE       1   		/* if photon not yet terminated */
#define DEAD        0    		/* if photon is to be terminated */
#define THRESHOLD   0.01		/* used in roulette */
#define CHANCE      0.1  		/* used in roulette */
#define Boolean     char
#define SQR(x)		(x*x)
#define SIGN(x)     ((x)>=0 ? 1:-1)
#define RandomNum   (double) RandomGen(1, 0, NULL) /* Calls for a random number. */
#define COS90D      1.0E-6      /* If cos(theta) <= COS90D, theta >= PI/2 - 1e-6 rad. */
#define ONE_MINUS_COSZERO 1.0E-12
/* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */


/* Propagation parameters */
double	x, y, z;        /* photon position */
double	ux, uy, uz;     /* photon trajectory as cosines */
double  uxx, uyy, uzz;	/* temporary values used during SPIN */
double	s;              /* step sizes. s = -log(RND)/mus [cm] */
double  sleft;          /* dimensionless */
double	costheta;       /* cos(theta) */
double  sintheta;       /* sin(theta) */
double	cospsi;         /* cos(psi) */
double  sinpsi;         /* sin(psi) */
double	psi;            /* azimuthal angle */
long	i_photon;       /* current photon */
double	W;              /* photon weight */
double	absorb;         /* weighted deposited in a step due to absorption */
short   photon_status;  /* flag = ALIVE=1 or DEAD=0 */
Boolean sv;             /* Are they in the same voxel? */

/* other variables */
double	mua;            /* absorption coefficient [cm^-1] */
double	mus;            /* scattering coefficient [cm^-1] */
double	g;              /* anisotropy [-] */
double	nr;              /* refractive index [-] */
double	mr;              /* mirror reflection chances [-] */
double	n1;              /* refractive index previous layer [-] */
double	n2;              /* refractive index next layer [-] */
double	Nphotons;       /* number of photons in simulation */

/* launch parameters */
int		mcflag, launchflag, boundaryflag;
float	xfocus, yfocus, zfocus;
float	ux0, uy0, uz0;
float	radius;
float	waist;

/* dummy variables */
double  rnd;            /* assigned random value 0-1 */
double	r, phi;			/* dummy values */
long	i,tempi,j,NN,Nyx,m;         /* dummy indices */
double	tempx, tempy, tempz; /* temporary variables, used during photon step. */
int 	ix, iy, iz;     /* Added. Used to track photons */
double 	temp;           /* dummy variable */
int     bflag;  // boundary flag; KE: 1 = photon inside volume ; 0 = outside
int 	surfflag; /* surface flag: 0 = photon inside tissue, 1 = escaped outside tissue */

/* mcxyz bin variables */
float	dx, dy, dz;     /* bin size [cm] */
int		Nx, Ny, Nz, Nt; /* # of bins */
float	xs, ys, zs;		/* launch position */
float 	zsurf, Rd;

/* time */
float	time_min;               // Requested time duration of computation.
time_t	now;
double	start_time, finish_time, temp_time; /* for clock() */
float   tsave; //RMT
float   tsaves; //RMT
float   t0; //RMT
float   t1s; //RMT
float   t1; //RMT
float   t2s; //RMT
float   t2; //RMT
float   t3s; //RMT
float   t3; //RMT

/* tissue parameters */
char	tissuename[50][32];
float 	muav[Ntiss];    // muav[0:Ntiss-1], absorption coefficient of ith tissue type
float 	musv[Ntiss];    // scattering coeff.
float 	gv[Ntiss];      // anisotropy of scattering
float 	nrv[Ntiss];      // refractive index
float 	mrv[Ntiss];      // mirror reflection chances

/**** KE start: Declaration of variables ****/
int face_dir; // exited voxel direction of the photon
int reflected; // check if the photon has been reflected
int det_num; // photon not detected yet/
double first_bias_done ; // photon not biased back - scattered yet
double cont_exist; // no split generated yet // check if a continuing photon packet exists
double L_current; // photon 's initial likelihood
double s_total; // photon 's initial path length
double z_max; // photon 's initial depth reached
int Ndetectors; // Number of source/detector pairs to represent conducting A-scans
int Pick_det; // index of randomly picked source/detector pair
double detx, dety, detz; // position of detector
double det_radius; // radius of detector
double cos_accept; // acceptance angle
//double a; //RMT no longuer useful
double costheta_S;
double costheta_B;
double sintheta_B;
double vx, vy, vz;
double upx, upy, upz;
double L_cont;
double i_cont;
double W_cont;
double x_cont, y_cont, z_cont;
double ux_cont, uy_cont, uz_cont;
double s_total_cont;
double z_max_cont;
double a_coef; // RMT
double p; // RMT
double det_z;
double f_HG, f_B;
long c_photon; // count collected photons
int *DetID;
float *DetW, *DetL, *DetS, *DetS2, *DetZ, *DetE; // RMT: I added DetE as a debugging variable for DetS Removed DetS since included in DetS
/**** KE end : Declaration of variables ****/


int main(int argc, const char * argv[])
{
    printf("argc = %d\n",argc);
    if (argc==0)
    {
		printf("which will load the files name_H.mci and name_T.bin\n");
		printf("and run the Monte Carlo program.\n");
		printf("Yields  name_F.bin, which holds the fluence rate distribution.\n");
        return 0;
    }

	/* Input/Output */
	//KE: char 	dirname[STRLEN];  // holds "../sims/"
	char   	myname[STRLEN];
	// Holds the user's choice of myname, used in input and output files.
	char	filename[STRLEN];     // temporary filename for writing output.
    FILE*	fid=NULL;             // file ID pointer
    char    buf[32];              // buffer for reading header.dat

    strcpy(myname, argv[1]);      // acquire name from argument of function call by user
    printf("name = %s\n",myname);

	/**** INPUT FILES *****/
    /* IMPORT myname_H.mci */
    strcpy(filename,myname);
    strcat(filename, "_H.mci");
	fid = fopen(filename,"r");
	fgets(buf, 32, fid);
	// run parameters
	sscanf(buf, "%lf", &Nphotons); // desired time duration of run [min] RMT now photons isntead
	fgets(buf, 32, fid);
	sscanf(buf, "%lf", &a_coef); // RMT chance of a foward photon doing a bias scattering.
	fgets(buf, 32, fid);
	sscanf(buf, "%lf", &p); // RMT chance of a foward photon doing a bias scattering.
	fgets(buf, 32, fid);
	sscanf(buf, "%d", &Ndetectors); // RMT chance of a foward photon doing a bias scattering.
	fgets(buf, 32, fid);
	sscanf(buf, "%d", &Nx);  // # of bins
	fgets(buf, 32,fid);
	sscanf(buf, "%d", &Ny);  // # of bins
	fgets(buf, 32,fid);
	sscanf(buf, "%d", &Nz);  // # of bins

	fgets(buf, 32,fid);
	sscanf(buf, "%f", &dx);	 // size of bins [cm]
	fgets(buf, 32,fid);
	sscanf(buf, "%f", &dy);	 // size of bins [cm]
	fgets(buf, 32,fid);
	sscanf(buf, "%f", &dz);	 // size of bins [cm]

	// launch parameters
	fgets(buf, 32,fid);
	sscanf(buf, "%d", &mcflag);  // mcflag, 0 = uniform, 1 = Gaussian, 2 = iso-pt
    // KE: 3 = rectangle
	fgets(buf, 32,fid);
	sscanf(buf, "%d", &launchflag);  // launchflag, 0 = ignore, 1 = manually set
    fgets(buf, 32,fid);
    sscanf(buf, "%d", &boundaryflag);  // 0 = no boundaries, 1 = escape at all boundaries, 2 = escape at surface only

	fgets(buf, 32,fid);
    sscanf(buf, "%f", &xs);  // initial launch point
	fgets(buf, 32,fid);
	sscanf(buf, "%f", &ys);  // initial launch point
	fgets(buf, 32,fid);
	sscanf(buf, "%f", &zs);  // initial launch point

	fgets(buf, 32,fid);
	sscanf(buf, "%f", &xfocus);  // xfocus
	fgets(buf, 32,fid);
	sscanf(buf, "%f", &yfocus);  // yfocus
	fgets(buf, 32,fid);
	sscanf(buf, "%f", &zfocus);  // zfocus

	fgets(buf, 32,fid);
	sscanf(buf, "%f", &ux0);  // ux trajectory
	fgets(buf, 32,fid);
	sscanf(buf, "%f", &uy0);  // uy trajectory
	fgets(buf, 32,fid);
	sscanf(buf, "%f", &uz0);  // uz trajectory

	fgets(buf, 32,fid);
	sscanf(buf, "%f", &radius);  // radius
	fgets(buf, 32,fid);
	sscanf(buf, "%f", &waist);  // waist

	fgets(buf, 32,fid);
	sscanf(buf, "%f", &zsurf);  // z_surface

	// tissue optical properties
	fgets(buf, 32,fid);
	sscanf(buf, "%d", &Nt);				// # of tissue types in tissue list
    printf("Nt = %d\n",Nt); // KE: check

	double s_total2[Nt]; // RMT : Create s_total here
	double s_total_cont2[Nt]; // RMT : Create s_total here

	for (i=1; i<=Nt; i++)
    {
        fgets(buf, 32, fid);
		sscanf(buf, "%f", &muav[i]);	// absorption coeff [cm^-1]
		fgets(buf, 32, fid);
		sscanf(buf, "%f", &musv[i]);	// scattering coeff [cm^-1]
		fgets(buf, 32, fid);
		sscanf(buf, "%f", &gv[i]);		// anisotropy of scatter [dimensionless]
		fgets(buf, 32, fid);
		sscanf(buf, "%f", &nrv[i]);		// refractive index
		fgets(buf, 32, fid);
		sscanf(buf, "%f", &mrv[i]);		// mirror reflection chances
	}
    fclose(fid);
	
	    printf("# of tissues available, Nt = %d\n",Nt);
    for (i=1; i<=Nt; i++)
    {
        printf("muav[%ld] = %0.4f [cm^-1]\n",i,muav[i]);
        printf("musv[%ld] = %0.4f [cm^-1]\n",i,musv[i]);
        printf("  gv[%ld] = %0.4f [--]\n",i,gv[i]);
		printf("  nrv[%ld] = %0.4f [--]\n",i,nrv[i]);
		printf("  mrv[%ld] = %0.4f [--]\n\n",i,mrv[i]);
    }
	
	
}