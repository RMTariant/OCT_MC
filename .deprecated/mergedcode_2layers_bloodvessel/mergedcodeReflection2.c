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

/* DECLARE FUNCTIONS */
double RandomGen(char Type, long Seed, long *Status);
/* Random number generator */
Boolean SameVoxel(double x1,double y1,double z1, double x2, double y2, double z2,
				double dx,double dy,double dz);
/* Asks,"In the same voxel?" */
double max2(double a, double b);
double min2(double a, double b);
double min3(double a, double b, double c);
double myErfInv2(double rndE);
double ReflectRefraction(double n1, double n2, double* x, double* y, double* z, double* ux, double* uy,
					double* uz, int* face_dir, int* reflected);
double MirrorReflection(double n1, double n2, double mr, double* x, double* y, double* z, double* ux, 
					double* uy, double* uz, int* face_dir, int* reflected);
double FindVoxelFace(double x1,double y1,double z1, double x2, double y2, double z2,
					double dx,double dy,double dz, double ux, double uy, double uz);
double FindVoxelFace2(double x1, double y1, double z1, int* det_num, int Pick_det,
        double detx, double det_radius, double det_z, double cos_accept,
        int Ndetectors, double dx, double dy, double dz, double ux, double uy, double uz, int* face_dir) ;

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
long	Nphotons;       /* number of photons in simulation */

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
	sscanf(buf, "%ld", &Nphotons); // desired time duration of run [min] RMT now photons isntead
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

    printf("Number of photons = %ld \n",Nphotons);
    printf("Nx = %d, dx = %0.4f [cm]\n",Nx,dx);
    printf("Ny = %d, dy = %0.4f [cm]\n",Ny,dy);
    printf("Nz = %d, dz = %0.4f [cm]\n",Nz,dz);

    printf("xs = %0.4f [cm]\n",xs);
    printf("ys = %0.4f [cm]\n",ys);
    printf("zs = %0.4f [cm]\n",zs);
    printf("mcflag = %d\n",mcflag);
    if (mcflag==0) printf("launching uniform flat-field beam\n");
    if (mcflag==1) printf("launching Gaissian beam\n");
    if (mcflag==2) printf("launching isotropic point source\n");
    if (mcflag==3) printf("launching square source\n");
    printf("xfocus = %0.4f [cm]\n",xfocus);
    printf("yfocus = %0.4f [cm]\n",yfocus);
    printf("zfocus = %0.2e [cm]\n",zfocus);
	if (launchflag==1)
	{
		printf("Launchflag ON, so launch the following:\n");
		printf("ux0 = %0.4f [cm]\n",ux0);
		printf("uy0 = %0.4f [cm]\n",uy0);
		printf("uz0 = %0.4f [cm]\n",uz0);
	}
	else
	{
		printf("Launchflag OFF, so program calculates launch angles.\n");
		printf("radius = %0.4f [cm]\n",radius);
		printf("waist  = %0.4f [cm]\n",waist);
	}
    printf("zsurf = %0.4f [cm]\n",zsurf);

    if (boundaryflag==0)
		printf("boundaryflag = 0, so no boundaries.\n");
    else if (boundaryflag==1)
		printf("boundaryflag = 1, so escape at all boundaries.\n");
	else if (boundaryflag==2)
		printf("boundaryflag = 2, so escape at surface only.\n");
	else
	{
        printf("improper boundaryflag. quit.\n");
        return 0;
    }
    printf("# of tissues available, Nt = %d\n",Nt);
    for (i=1; i<=Nt; i++)
    {
        printf("muav[%ld] = %0.4f [cm^-1]\n",i,muav[i]);
        printf("musv[%ld] = %0.4f [cm^-1]\n",i,musv[i]);
        printf("  gv[%ld] = %0.4f [--]\n",i,gv[i]);
		printf("  nrv[%ld] = %0.4f [--]\n",i,nrv[i]);
		printf("  mrv[%ld] = %0.4f [--]\n\n",i,mrv[i]);
    }

    /* IMPORT BINARY TISSUE FILE */
	char 	*v=NULL;
	float 	*F=NULL;
	float	*R=NULL;
	int 	type;
	NN = Nx*Ny*Nz;
	Nyx = Nx*Ny;
	v  = ( char *)malloc(NN*sizeof(char));  /* tissue structure */
	F  = (float *)malloc(NN*sizeof(float));	/* relative fluence rate [W/cm^2/W.delivered] */
	R  = (float *)malloc(Nyx*sizeof(float));	/* escaping flux [W/cm^2/W.delivered] */

    DetID  = malloc(sizeof(int));	// KE: photon ID at det_num
    DetS = malloc(sizeof(float)); // KE: photon path length
	DetS2 = malloc(sizeof(float)); // KE: photon path length
	DetE = malloc(sizeof(float)); // RMT: Debugging variable. Plays no role.
    DetW  = malloc(sizeof(float));	// KE: photon weight
    DetL  = malloc(sizeof(float));	// KE: likelihood ratio
    DetZ  = malloc(sizeof(float));	// KE: photons reached depth
	// read binary file
    strcpy(filename,myname);
    strcat(filename, "_T.bin");
    fid = fopen(filename, "rb");
    fread(v, sizeof(char), NN, fid);
    fclose(fid);

    // Show tissue on screen, along central z-axis, by listing tissue type #'s.
    iy = Ny/2;
    ix = Nx/2;
    printf("central axial profile of tissue types:\n");
    for (iz=0; iz<Nz; iz++)
    {
        i = (long)(iz*Ny*Nx + ix*Ny + iy);
        printf("%d",v[i]);
    }
    printf("\n\n");

	/**************************
	 * ============================ MAJOR CYCLE ========================
	 **********/
	start_time = clock();
	now = time(NULL);
	printf("%s\n", ctime(&now));

	/**** INITIALIZATIONS
	 *****/
	RandomGen(0, -(int)time(NULL)%(1<<15), NULL);
	/* initiate with seed = 1, or any long integer. */
	for(j=0; j<NN;j++) 	F[j] = 0.0; // ensure F[] starts empty.
	for(j=0; j<Nyx;j++) R[j] = 0.0;  // RMT ensure R[] starts empty.
	Rd = 0.0;

	/**** RUN
	 Launch N photons, initializing each one before progation.
	 *****/
	printf("------------- Begin Monte Carlo -------------\n");
    printf("%s\n",myname);
    printf("requesting unknown time\n");
	//Nphotons = 200; // will be updated to achieve desired run time, time_min.
	i_photon = 0;
    c_photon = 0;
    //a = 0.925; //KE: Lima et al 2012
	do {
		//t0 = clock(); //RMT
        // KE: while (i_photon < Nphotons) RMT: Main loop simulating all photons.
		/**** LAUNCH: Initialize photon position and trajectory *****/
		i_photon += 1;				/* increment photon count */
		W = 1.0;                    /* set photon weight to one */
        //printf("W = %f\n",W); // KE: check W
		photon_status = ALIVE;      /* Launch an ALIVE photon */

        /*** KE start: this part is from the A.3 of Zhao's thesis ***/
        det_num = -1; /* photon not detected yet */
        first_bias_done = 0; /* photon not biased back - scattered yet */
        cont_exist = 0; /* no split generated yet */
        L_current = 1; /* photon 's initial likelihood */
        s_total = 0; /* photon 's initial path length */
		for (m=0; m<Nt; m++)
		{
			s_total2[m] = 0;
		}
        z_max = 0; /* photon 's initial depth reached */
        //Ndetectors = 512; // KE: number of detectors RMT: Changed to one for debugging
        det_radius = 0.001; //KE: Zhao's thesis chapter 3.3.4 10micro m
        cos_accept = cos(5); //KE: Zhao's thesis chapter 3.3.4

        //a = (double) RandomGen(1, 0, NULL); //KE


        /* pick the fiber that the current photon packet is biased towards to: [1, Ndetectors ] */
        // KE: array of source/detector pairs
        while ((rnd = RandomGen(1, 0, NULL)) >= 1.0) ; // avoids rnd = 1
        Pick_det = floor(rnd * Ndetectors) + 1; // KE: randomly pick one source/detector pair
        //printf("Pick_det = %d\n",Pick_det); //KE: check Pick_det

        /* Set trajectory to mimic A- scan */
        // KE: each A-scan launches and collects photon packets at a different lateral position.
        // KE: That is the only difference among A-scans
        if (Ndetectors == 1)
        {
            detx = 0;
        }
        else
        {
            detx = 2 * radius * (Pick_det - 1) / (Ndetectors - 1) - radius;
        }
        //printf("detx = %f\n",detx); //KE: check detx

        /*** KE end : this part is from the A.3 of Zhao's thesis ***/

		// Print out message about progress.
		if ((i_photon>200) & (fmod(i_photon, (int)(Nphotons/20))  == 0))
		{
            printf("%ld%% done\n", i_photon/Nphotons*100); //RMT
			//temp = i_photon/Nphotons*100;
            //printf("%0.1f%% \t\tfmod = %0.3f\n", temp,fmod(temp, 10.0));
            //if ((temp<10) | (temp>90))
            //{
            //    printf("%0.0f%% done\n", i_photon/Nphotons*100);
            //}
            //else if(fmod(temp, 10.0)>9) //RMT
            //    printf("%0.0f%% done\n", i_photon/Nphotons*100);
        }

		// At 1000th photon, update Nphotons to achieve desired runtime (time_min)
		if (i_photon==1)
			temp_time = clock();
		if (i_photon==200)
		{
			finish_time = clock();
			//Nphotons = (long)( time_min*60*999*CLOCKS_PER_SEC/(finish_time-temp_time) );
			printf("Nphotons = %ld for unkown time\n",Nphotons)
			;
		}
		/**** SET SOURCE
		 * Launch collimated beam at x,y center.
		 ****/

		/****************************/
		/* Initial position. */

		/* trajectory */ //
		if (launchflag==1) // manually set launch
		{
            x	= xs + detx;
			y	= ys;
			z	= zs;
			ux	= 0;
			uy	= 0;
			uz	= 1;
		}
		else // use mcflag
		{
			if (mcflag==0)
			{ // uniform beam
                // set launch point and width of beam
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				r		= radius*sqrt(rnd); // radius of beam at launch point
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				phi		= rnd*2.0*PI;
				x		= xs + r*cos(phi);
				y		= ys + r*sin(phi);
				z		= zs;
				// set trajectory toward focus
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				r		= waist*sqrt(rnd); // radius of beam at focus
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				phi		= rnd*2.0*PI;
				xfocus	= r*cos(phi);
				yfocus	= r*sin(phi);
				temp	= sqrt((x - xfocus)*(x - xfocus) + (y - yfocus)*(y - yfocus)
					            + zfocus*zfocus);
				ux		= -(x - xfocus)/temp;
				uy		= -(y - yfocus)/temp;
				uz		= sqrt(1 - ux*ux - uy*uy);
			}
			else if (mcflag==1)
			{ // gaussian beam with a 1/e2 waist
				ux	= 0;
				uy	= 0;
				uz	= 1;
				z	= zs;
				y   = ys;
				rnd = myErfInv2(2*RandomGen(1,0,NULL)-1);
				x   = xs + detx + (waist*rnd/2/sqrt(2));

			}
			else if (mcflag==2)
			{ // isotropic pt source
				costheta = 1.0 - 2.0*RandomGen(1,0,NULL);
				sintheta = sqrt(1.0 - costheta*costheta);
				psi = 2.0*PI*RandomGen(1,0,NULL);
				cospsi = cos(psi);
				if (psi < PI)
					sinpsi = sqrt(1.0 - cospsi*cospsi);
				else
					sinpsi = -sqrt(1.0 - cospsi*cospsi);
				x = xs;
				y = ys;
				z = zs;
				ux = sintheta*cospsi;
				uy = sintheta*sinpsi;
				uz = costheta;
			}
			else if (mcflag==3)
			{ // rectangular source collimated
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				x = radius*(rnd*2-1); // use radius to specify x-halfwidth of rectangle
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				y = radius*(rnd*2-1); // use radius to specify y-halfwidth of rectangle
				z = zs;
				ux = 0.0;
				uy = 0.0;
				uz = 1.0; // collimated beam
			}
		} // end  use mcflag
		/****************************/

		/* Get tissue voxel properties of launchpoint.
		 * If photon beyond outer edge of defined voxels,
		 * the tissue equals properties of outermost voxels.
		 * Therefore, set outermost voxels to infinite background value.
		 */
		ix = (int)(Nx/2 + x/dx);
		iy = (int)(Ny/2 + y/dy);
		iz = (int)(z/dz);
		if (ix>=Nx) ix=Nx-1;
		if (iy>=Ny) iy=Ny-1;
		if (iz>=Nz) iz=Nz-1;
		if (ix<0)   ix=0;
		if (iy<0)   iy=0;
		if (iz<0)   iz=0;
		/* Get the tissue type of located voxel */
		i		= (long)(iz*Ny*Nx + ix*Ny + iy);
		type	= v[i];
		mua 	= muav[type];
		mus 	= musv[type];
		g 		= gv[type];
		nr 		= nrv[type];
		mr 		= mrv[type];

        bflag = 1;
        // initialize as 1 = inside volume, but later check as photon propagates.
        surfflag = 1; // initially inside tissue
        // NOTE: must launch photons at tissue surface, or surfflag will go to 0.
        det_z = z; //KE
		
		//t1 = clock()-t0; //RMT
		//t1s += t1; //RMT

		/* HOP_DROP_SPIN_CHECK
		 Propagate one photon until it dies as determined by ROULETTE.
		 *******/
		do {
			/**** HOP
			 Take step to new position
			 s = dimensionless stepsize
			 x, uy, uz are cosines of current photon trajectory
			 *****/
			while ((rnd = RandomNum) <= 0.0);   /* yields 0 < rnd <= 1 */
			sleft	= -log(rnd);				/* dimensionless step */
			//KE CNT += 1;
			// printf("sleft = %f\n",sleft); // KE: check sleft


			if (photon_status == DEAD) //
            { // load the continuing photon and update the flags

                 x = x_cont;
                 y = y_cont;
                 z = z_cont;
                 ux = ux_cont;
                 uy = uy_cont;
                 uz = uz_cont;
                 i = i_cont;
				 s_total = s_total_cont;
				 for (m=0; m<Nt; m++)
				{
					s_total2[m] = s_total_cont2[m];
				}
                 
                 z_max = z_max_cont;
                 type = v[i];
                 mua = muav[type];
                 mus = musv[type];
                 g = gv[type];
				 nr = nrv[type];
				 mr = mrv[type];
                 W = W_cont;
                 L_current = L_cont;
                 cont_exist = 0;
                 photon_status = ALIVE;
                 first_bias_done = 0;
                 det_num = -1;
            }


			do{  // while sleft>0
			// RMT: in this loop, we are looking at the different medium the photon is
			//      going through. Where it is absorbed, detected or escape the simulation
                s     = sleft/mus;				/* Step size [cm].*/
                // printf("s = %f\n",s); // KE: check mus

				tempx = x + s*ux;				/* Update positions. [cm] */
				tempy = y + s*uy;
				tempz = z + s*uz;

				sv = SameVoxel(x,y,z, tempx, tempy, tempz, dx,dy,dz);
				if (sv) /* photon in same voxel */
				{
                    //printf("photon in same voxel\n");
					x=tempx;					/* Update positions. */
					y=tempy;
					z=tempz;

					/**** DROP
					 Drop photon weight (W) into local bin.
					 *****/
                    absorb = W*(1 - exp(-mua*s));
                    /* photon weight absorbed at this step */
                    W -= absorb;
                    //printf("W = %f\n",W); // KE: check W
                    /* decrement WEIGHT by amount absorbed */
					// If photon within volume of heterogeneity, deposit energy in F[].
					// Normalize F[] later, when save output.
                    if (bflag) F[i] += absorb*L_current;
                    // only save data if blag==1, i.e., photon inside simulation cube

					/* Update sleft */
					sleft = 0;		/* dimensionless step remaining */

				    /* Update total path length */ // RMT
					s_total += s;
				    s_total2[type-1] += s; // RMT Update the total distance here.
				}
				else /* photon has crossed voxel boundary */
				{
					/* wakka */
					
					
					/* step to voxel face + "littlest step" so just inside new voxel. */
                    s = ls + FindVoxelFace2(x, y, z, &det_num, Pick_det, detx, det_radius, det_z, cos_accept, Ndetectors, dx, dy, dz, ux, uy, uz, &face_dir);
                    //s_total += s; // RMT Update the total distance here. Not suppose to be here

					/*** DROP: Drop photon weight (W) into local bin  ***/
					absorb = W*(1-exp(-mua*s));  /* photon weight absorbed at this step */
					W -= absorb;                 /* decrement WEIGHT by amount absorbed */

					// If photon within volume of heterogeneity, deposit energy in F[].
					// Normalize F[] later, when save output.
                    if (bflag) F[i] += absorb*L_current;

                    if (det_num != -1)
                    { /* check if the photon is detected . */
                        // KE: det_num changes in FindVoxelFace2 function when photon gets detected

                      /* Update total path length */
                      s_total += s;
					  s_total2[type-1] += s;

							//tsave = clock();
                         /* Save properties of interest */
                         if (L_current > 0 &&  det_num == Pick_det)
                         { // avoid NAN and zero likelihood, and avoid cross - detection
                             // Def: float *DetW, *DetL, *DetS, *DetZ;
                             // DetS  = malloc(sizeof(float));
                             DetS = realloc(DetS,(c_photon+2)* sizeof(float)); //RMT
							 DetS2 = realloc(DetS2,((c_photon+2)*Nt)* sizeof(float));
							 DetE = realloc(DetE,(c_photon+2)* sizeof(float)); // RMT
                             DetS[c_photon]=s_total; //RMT
							 for (m=0; m<Nt; m++)
							 {
								 DetS2[Nt*c_photon+m] = s_total2[m];
							 }
							 DetE[c_photon]=a_coef; // RMT, redone
                             DetID = realloc(DetID,(c_photon+2)* sizeof(int));
                             DetID[c_photon] = det_num;
                             DetW = realloc(DetW,(c_photon+2)* sizeof(float));
                             DetW[c_photon] = W;
                             DetL= realloc(DetL,(c_photon+2)* sizeof(float));
                             DetL[c_photon] = L_current;
                             DetZ = realloc(DetZ,(c_photon+2)* sizeof(float));
                             DetZ[c_photon] = z_max;
                             /* increment collected photon count */
                             c_photon += 1;
                         }
							//tsave -= clock(); //RMT
							//tsaves += tsave; //RMT
                         // if( c_photon ==1) { printf (" OK at 590;\ n") ;}
                         photon_status = DEAD; // RMT This might need to be "dead"
                         sleft = 0;
                    }
                    else // RMT photon not detected
                    {
                        /* Update sleft */
                        sleft -= s*mus;  /* dimensionless step remaining */
                        if (sleft<=ls) sleft = 0;

                        /* Update positions. */
                        x += s*ux;
                        y += s*uy;
                        z += s*uz;

						/* Update total path length */ // RMT
						s_total += s; //RM
						s_total2[type-1] += s; //RMT

                        // pointers to voxel containing optical properties
                        ix = (int)(Nx/2 + x/dx);
                        iy = (int)(Ny/2 + y/dy);
                        iz = (int)(z/dz);
                        if (ix>Nx) ix=Nx;
                        if (iy>Ny) iy=Ny;
                        if (ix<0)  ix=0;
                        if (iy<0)  iy=0;


                        //*** ESCAPE or not
                        if((surfflag==1) & (z<=zsurf)) // escape at surface
                        {
                            Rd += W;
                            i = (long)(Nx*ix + iy);
                            R[i] += W;
                            surfflag = 0;  // disable repeated assignments to Rd, R[i]
                        }
                        if (z<0) // escape cube
                        {
                            photon_status = DEAD; // RMT: as in really dead this time
                            sleft = 0;
                        }
                        else // No escape
                        {
                            bflag = 1;
                            // Boundary flag. Initialize as 1 = inside volume, then check.
                            if (boundaryflag==0) // Infinite medium
                            {
                            // Check if photon has wandered outside volume.
                            // If so, set tissue type to boundary value, but let photon wander
                            // Set blag to zero, so DROP does not deposit energy.
                                if (iz>=Nz) {iz=Nz-1; bflag = 0;}
                                if (ix>=Nx) {ix=Nx-1; bflag = 0;}
                                if (iy>=Ny) {iy=Ny-1; bflag = 0;}
                                if (iz<0)   {iz=0;    bflag = 0;}
                                if (ix<0)   {ix=0;    bflag = 0;}
                                if (iy<0)   {iy=0;    bflag = 0;}
                            }
                            else if (boundaryflag==1) // Escape at boundaries
                            {
                                if (iz>=Nz) {iz=Nz-1; photon_status = DEAD; sleft=0;}
                                if (ix>=Nx) {ix=Nx-1; photon_status = DEAD; sleft=0;}
                                if (iy>=Ny) {iy=Ny-1; photon_status = DEAD; sleft=0;}
                                if (iz<0)   {iz=0;    photon_status = DEAD; sleft=0;}
                                if (ix<0)   {ix=0;    photon_status = DEAD; sleft=0;}
                                if (iy<0)   {iy=0;    photon_status = DEAD; sleft=0;}
                            }
                            else if (boundaryflag==2)
                            { // Escape at top surface, no x,y bottom z boundaries
                                if (iz>=Nz) {iz=Nz-1; bflag = 0;}
                                if (ix>=Nx) {ix=Nx-1; bflag = 0;}
                                if (iy>=Ny) {iy=Ny-1; bflag = 0;}
                                if (iz<0)   {iz=0;    photon_status = DEAD; sleft=0;}
                                if (ix<0)   {ix=0;    bflag = 0;}
                                if (iy<0)   {iy=0;    bflag = 0;}
                            }
                            // update pointer to tissue type
							tempi = i;
                            i    = (long)(iz*Ny*Nx + ix*Ny + iy);
							n1 = nr;
							type = v[i];
							n2 = nrv[type];
							mr = mrv[type];
							reflected = 0;
							// RMT check if reflection
							if (mr != 0 && face_dir == 3)
							{
								// RMT add mirror reflection
								MirrorReflection(n1, n2, mr, &x, &y, &z, &ux, &uy, &uz, &face_dir, &reflected);
							}
							else
							{
								// RMT check if change in index of refraction and apply reflection/refraction
								if(n1 != n2) 
								{
									ReflectRefraction(n1, n2, &x, &y, &z, &ux, &uy, &uz, &face_dir, &reflected);
								}
							}
							// Update the optical properties depending if in next or same voxel
							if(reflected == 1)
							{
								i = tempi;
							}
							type = v[i];
							mua  = muav[type];
							mus  = musv[type];
							g    = gv[type];
							nr	 = nrv[type];
							mr	 = mrv[type];

                        }
                    }
				} //(sv) /* same voxel */

			} while(sleft>0); //do...while

			//t2 = clock() - t1; //RMT
			//t2s += t2; //RMT

            /***
            * KE start: this part is from the A.4 of Zhao's thesis
            ***/
            /**** SPIN AND SPLIT
            * The Spin process is to scatter photon into new
            trajectory defined by theta and psi. Theta is specied by cos(theta) , which is determined based
            on the Henyey-Greenstein scattering function, and then convert theta and psi into
            cosines ux, uy, uz. Split follows exactly the procedure described in Section 4.3,
            where we apply biased backward-scatterings and biased forward-scattering, as well
            as unbiased scatterings. Once therst biased backward-scattering takes place, we
            split the photon packet into two if the likelihood ratio of the biased back-scattering
            is less than 1. We save the information of the continuing photon and continue to
            track the current photon, by applying biased and unbiased forward-scatterings*/

            if (photon_status == ALIVE)
            {
                 /* check whether the first biased back - scattering has been applied : 0 = not ...
                  applied , 1 = applied */


                 if (first_bias_done == 0)
                { /* apply the first biased scattering */
                     /* Sample for costheta_B */
                     rnd = RandomNum;
                     //a = RandomNum; // KE: added RMT WRONG! IT IS A PARAMETER!
                     double temp = 1 / sqrt(1 + a_coef*a_coef) + rnd * (1 / (1 - a_coef) - 1 / sqrt(1 + a_coef*a_coef));
                     costheta_B = 0.5 / a_coef*(a_coef*a_coef + 1 - 1 / (temp*temp));
                     sintheta_B = sqrt(1.0 - costheta_B * costheta_B);
                     /* Sample psi . */
                     psi = 2.0 * PI * RandomNum;
                     cospsi = cos(psi);
                     if (psi < PI)
                         sinpsi = sqrt(1.0 - cospsi * cospsi); /* sqrt () is faster than sin (). */
                     else
                         sinpsi = -sqrt(1.0 - cospsi * cospsi);
                     /* Compute the unit vector v towards the actual position of the detector , ...
                      where detx is chosen uniformly along the centers of the collecting fiber ...
                      array . */
                     // KE: biased direction; direction of the actual position  of the collecting optics
                     if (Ndetectors == 1)
                     {
                         detx = 0;
                     }
                     else
                     {
                         detx = 2 * radius * (Pick_det - 1) / (Ndetectors - 1) - radius;
                     }
                     dety = 0;
                     detz = det_z;
                     temp = sqrt((x - detx) * (x - detx) + (y - dety) * (y - dety) + (z - detz) * (z - detz));
                     vx = -(x - detx) / temp;
                     vy = -(y - dety) / temp;
                     vz = -(z - detz) / temp;

                     /* New trajectory u' = (upx , upy , upz) */
                     // KE: equal to equation 3.21 in the mcmcl manualn of Steves Jacques
                     if (1 - fabs(vz) <= ONE_MINUS_COSZERO)
                     { /* close to perpendicular . */
                         upx = sintheta_B * cospsi;
                         upy = sintheta_B * sinpsi;
                         upz = costheta_B * SIGN(vz); /* SIGN () is faster than division . */
                     }
                     else
                     {                      /* usually use this option */
                         // KE: equal to equation 3.22 in the mcmcl manualn of Steves Jacques
                         temp = sqrt(1.0 - vz * vz);
                         upx = sintheta_B * (vx * vz * cospsi - vy * sinpsi) / temp + vx * costheta_B;
                         upy = sintheta_B * (vy * vz * cospsi + vx * sinpsi) / temp + vy * costheta_B;
                         upz = -sintheta_B * cospsi * temp + vz * costheta_B;
                     }
                     /* Compute the likelihood ratio for this particular biased ...
                      back - scattering */
                     // KE: henyey greenstein probability; equal to equation 4.3 in Zhao's thesis

                     costheta_S = upx * ux + upy * uy + upz * uz;
                     temp = (1 + a_coef * a_coef - 2 * a_coef * costheta_B) / (1 + g * g - 2 * g * costheta_S);
                     double L_temp = (1 - g * g) / (2 * a_coef * (1 - a_coef)) * (1 - (1 - a_coef) / sqrt(1 + a_coef * a_coef)) * sqrt(temp * temp * temp);

                     /* Check do we have a continuing photon packet or not ? */
                     // KE: this part is explained in section 4.3.3 in Zhao's thesis
                     if (L_temp < (1 - ls))
                     { // yes , do the unbiased spin and save the trajectory for the continuing photon packet
                         L_cont = L_current * (1 - L_temp);
                         i_cont = i;
                         /* the unbiased spin */
                         // KE: equal to equation 3.19 in the manual
                         /* Sample for costheta */
                         rnd = RandomNum;
                         if (g == 0.0)
                         costheta = 2.0 * rnd - 1.0;
                         else
                         {
                           double temp = (1.0 - g * g) / (1.0 - g + 2 * g * rnd);
                           costheta = (1.0 + g * g - temp * temp) / (2.0 * g);
                         }
                         sintheta = sqrt(1.0 - costheta * costheta); /* sqrt () is faster than sin (). */
                         /* Sample psi . */
                         // KE: equal to equation 3.22 in the manual RMT: From the original code
                         psi = 2.0 * PI * RandomNum;
                         cospsi = cos(psi);
                         if (psi < PI)
                                sinpsi = sqrt(1.0 - cospsi * cospsi); /* sqrt () is faster than sin (). */
                         else
                                sinpsi = -sqrt(1.0 - cospsi * cospsi);
                         /* New trajectory . */
                         // KE: equal to equation 3.22 in the manual RMT: From the original code
                         if (1 - fabs(uz) <= ONE_MINUS_COSZERO)
                         { /* close to perpendicular . */
                           uxx = sintheta * cospsi;
                           uyy = sintheta * sinpsi;
                           uzz = costheta * SIGN(uz); /* SIGN () is faster than division . */
                         }
                         else
                         { /* usually use this option */
                           // KE: equal to equation 3.21 in the manual RMT: From the original code
                           temp = sqrt(1.0 - uz * uz);
                           uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
                           uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
                           uzz = -sintheta * cospsi * temp + uz * costheta;
                         }
                         // KE: update of continuing photon direction
                         ux_cont = uxx;
                         uy_cont = uyy;
                         uz_cont = uzz;
                         // KE: save information of the contunuing photon
                         x_cont = x;
                         y_cont = y;
                         z_cont = z;
                         W_cont = W;
                         s_total_cont = s_total;
						 for (m=0; m<Nt; m++)
						 {
						 	s_total_cont2[m] = s_total2[m];
						 }
                         z_max_cont = z_max;
                         L_current *= L_temp;
                         cont_exist = 1;
                    }
                    else
                     { // no continuing photon packet
                       // KE: no unbiased spin after first biased scattering
                         L_current *= L_temp;
                         cont_exist = 0;
                     }
                     /* Update trajectory */
                     // KE: keep simulating the back scattered photon packet
                     ux = upx;
                     uy = upy;
                     uz = upz;
                     first_bias_done = 1;
                }
          // KE: chapter 4.3.3 in Zhao's thesis
				else
				{/* first biased back - scattering already done , apply additional biased ...
              forward - scattering */
					if (RandomNum <= p)
					{ // apply biased forward - scattering
                  /* Sample for costheta_B */
                  rnd = RandomNum;
                  double temp = 1 / sqrt(1 + a_coef * a_coef) + rnd * (1 / (1 - a_coef) - 1 / sqrt(1 + a_coef * a_coef));
                  costheta_B = 0.5 / a_coef * (a_coef * a_coef + 1 - 1 / (temp * temp));
                  sintheta_B = sqrt(1.0 - costheta_B * costheta_B);
                  /* Sample psi . */
                  psi = 2.0 * PI * RandomNum;
                  cospsi = cos(psi);
                  if (psi < PI)
                      sinpsi = sqrt(1.0 - cospsi * cospsi); /* sqrt () is faster than sin (). */
                  else
                      sinpsi = -sqrt(1.0 - cospsi * cospsi);
                  /* Compute the unit vector v towards the actual position of the ...
                  detector , where detx is chosen uniformly along the centers of the ...
                  collecting fiber array . */
                  if (Ndetectors == 1)
                      detx = 0;
                  else
                      detx = 2 * radius * (Pick_det - 1) / (Ndetectors - 1) - radius; // RMT: detx ready calculated
                  dety = 0;
                  detz = det_z;
                  temp = sqrt((x - detx) * (x - detx) + (y - dety) * (y - dety) + (z - detz) * (z - detz));
                  vx = -(x - detx) / temp;
                  vy = -(y - dety) / temp;
                  vz = -(z - detz) / temp;
                  /* New trajectory u' = (upx , upy , upz) */
                  if (1 - fabs(vz) <= ONE_MINUS_COSZERO)
                  {/* close to perpendicular . */
                      upx = sintheta_B * cospsi;
                      upy = sintheta_B * sinpsi;
                      upz = costheta_B * SIGN(vz); /* SIGN () is faster than division . */
                  }
                  else
                  { /* usually use this option */
                      temp = sqrt(1.0 - vz * vz);
                      upx = sintheta_B * (vx * vz * cospsi - vy * sinpsi) / temp + vx * costheta_B;
                      upy = sintheta_B * (vy * vz * cospsi + vx * sinpsi) / temp + vy * costheta_B;
                      upz = -sintheta_B * cospsi * temp + vz * costheta_B;
                  }
                  // KE: equal to equation 4.3 in Zhao's thesis
                  /* Compute the likelihood ratio for this particular biased ...
                 forward - scattering */
                  costheta_S = upx * ux + upy * uy + upz * uz;
                  temp = 1 + g * g - 2 * g * costheta_S;
                  f_HG = (1 - g * g) * 0.5 / sqrt(temp * temp * temp);
                  temp = 1 + a_coef * a_coef - 2 * a_coef * costheta_B;
                  f_B = a_coef * (1 - a_coef) / ((sqrt(temp * temp * temp)) * (1 - (1 - a_coef) / sqrt(1 + a_coef * a_coef)));
                  double L_temp = f_HG / (p * f_B + (1 - p) * f_HG); // KE: equal to equation 4.4 in Zhao's thesis
                  L_current *= L_temp;
                  /* Update trajectory */
                  ux = upx;
                  uy = upy;
                  uz = upz;
					}
					else
					{// apply unbiased scattering
              /* Sample for costheta */
                  rnd = RandomNum;
                  if (g == 0.0)
                      costheta = 2.0 * rnd - 1.0;
                  else
                  {
                      double temp = (1.0 - g * g) / (1.0 - g + 2 * g * rnd);
                      costheta = (1.0 + g * g - temp * temp) / (2.0 * g);
                  }
                  sintheta = sqrt(1.0 - costheta * costheta); /* sqrt () is faster than sin (). */
                  /* Sample psi . */
                  psi = 2.0 * PI * RandomNum;
                  cospsi = cos(psi);
                  if (psi < PI)
                      sinpsi = sqrt(1.0 - cospsi * cospsi); /* sqrt () is faster than sin (). */
                  else
                      sinpsi = -sqrt(1.0 - cospsi * cospsi);
                  /* New trajectory . */
                  if (1 - fabs(uz) <= ONE_MINUS_COSZERO)
                  { /* close to perpendicular . */
                      uxx = sintheta * cospsi;
                      uyy = sintheta * sinpsi;
                      uzz = costheta * SIGN(uz); /* SIGN () is faster than division . */
                  }
                  else
                  { /* usually use this option */
                      temp = sqrt(1.0 - uz * uz);
                      uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
                      uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
                      uzz = -sintheta * cospsi * temp + uz * costheta;
                  }
                  /* Compute the unit vector v towards the actual position of the ...
                  detector , where detx is chosen uniformly along the centers of the ...
                  collecting fiber array . */
				  // RMT: I don't think calculating detx is necessary here. It was already done before.
                  if (Ndetectors == 1)
                      detx = 0;
                  else
                      detx = 2 * radius * (Pick_det - 1) / (Ndetectors - 1) - radius;
                  dety = 0;
                  detz = det_z;
                  temp = sqrt((x - detx) * (x - detx) + (y - dety) * (y - dety) + (z - detz) * (z - detz));
                  vx = -(x - detx) / temp;
                  vy = -(y - dety) / temp;
                  vz = -(z - detz) / temp;
                  /* Compute the likelihood ratio for this particular unbiased ...
                  forward - scattering */
                  costheta_S = costheta;
                  costheta_B = uxx * vx + uyy * vy + uzz * vz;
                  temp = 1 + g * g - 2 * g * costheta_S;
                  f_HG = (1 - g * g) * 0.5 / sqrt(temp * temp * temp);
                  temp = 1 + a_coef * a_coef - 2 * a_coef * costheta_B;
                  f_B = a_coef * (1 - a_coef) / ((sqrt(temp * temp * temp)) * (1 - (1 - a_coef) / sqrt(1 + a_coef * a_coef)));
                  double L_temp = f_HG / (p * f_B + (1 - p) * f_HG);
                  L_current *= L_temp;
                  /* Update trajectory */
                  ux = uxx;
                  uy = uyy;
                  uz = uzz;
					}
				}
				//t3 = clock()-t2;
				//t3s += t3;
			/***
			* KE end : this part is from the A.4 of Zhao's thesis
			***/
			/**** CHECK ROULETTE
	   If photon weight below THRESHOLD, then terminate photon using Roulette
	   technique. Photon has CHANCE probability of having its weight increased
	   by factor of 1/CHANCE, and 1-CHANCE probability of terminating.
	   *****/

				if (W < THRESHOLD)
				{
					if (RandomNum <= CHANCE)
						W /= CHANCE;
					else
						photon_status = DEAD;
				}
			}
		}	 while (photon_status == ALIVE || cont_exist == 1 );  /* end STEP_CHECK_HOP_SPIN */
        /* if ALIVE, continue propagating */
		/* If photon DEAD, then launch new photon. */

	} while (i_photon < Nphotons);  /* end RUN */
	//RMT: End of the monte carlo simulation. All photons were simulated.


	printf("collected photons = %ld\n",c_photon);

	printf("------------------------------------------------------\n");
	finish_time = clock();
	time_min = (double)(finish_time-start_time)/CLOCKS_PER_SEC/60;
	printf("Elapsed Time for %ld photons = %5.3f min\n",Nphotons,time_min);
	printf("%0.02e photons per minute\n", Nphotons/time_min);
	//printf("%0.2e saving files\n", tsaves/CLOCKS_PER_SEC/60);
	//printf("%0.2e set photon start\n", t1s/CLOCKS_PER_SEC/60);
	//printf("%0.2e propagate\n", t2s/CLOCKS_PER_SEC/60);
	//printf("%0.2e spin and split\n", t3s/CLOCKS_PER_SEC/60);
    /**** SAVE
     Convert data to relative fluence rate [cm^-2] and save.
     *****/

     // Save the binary file



    strcpy(filename,myname);
    strcat(filename,"_DetS.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(DetS, sizeof(float), c_photon, fid);
    fclose(fid);
	
    strcpy(filename,myname);
    strcat(filename,"_DetS2.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(DetS2, sizeof(float), c_photon*Nt, fid);
    fclose(fid);

	strcpy(filename,myname); // RMT extra saved data. Whole paragraph was added.
    strcat(filename,"_DetE.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(DetE, sizeof(float), c_photon, fid);
    fclose(fid);

     // Save the binary file
    strcpy(filename,myname);
    strcat(filename,"_DetW.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(DetW, sizeof(float), c_photon, fid);
    fclose(fid);

     // Save the binary file
    strcpy(filename,myname);
    strcat(filename,"_DetL.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(DetL, sizeof(float), c_photon, fid);
    fclose(fid);

     // Save the binary file
    strcpy(filename,myname);
    strcat(filename,"_DetZ.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(DetZ, sizeof(float), c_photon, fid);
    fclose(fid);

    // Save the binary file
    strcpy(filename,myname);
    strcat(filename,"_DetID.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(DetID, sizeof(int), c_photon, fid);
    fclose(fid);

    // Normalize deposition (A) to yield fluence rate (F).
    temp = dx*dy*dz*Nphotons;
    for (i=0; i<NN;i++) F[i] /= (temp*muav[v[i]]);
    // Save the binary file
    strcpy(filename,myname);
    strcat(filename,"_F.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(F, sizeof(float), NN, fid);
    fclose(fid);

    /* save reflectance */ // RMT: Not fully implemented. I don't think it is useful to implement it.
    temp = dx*dy*Nphotons;
	for (i=0; i<Nyx; i++) R[i] /= temp;
	strcpy(filename,myname);
	strcat(filename,"_Ryx.bin");
	printf("saving %s\n",filename);
	fid = fopen(filename, "wb");   /* 2D voxel output */
	fwrite(R, sizeof(float), Nyx, fid);
	fclose(fid);

	printf("WRd = %0.0f\n",Rd);
	printf("Nphotons = %ld\n",Nphotons);
	Rd /= Nphotons;
	printf("Rd = %0.3e\n",Rd);
	strcpy(filename,myname);
	strcat(filename,"_Rd.dat");
	printf("saving %s\n",filename);
	fid = fopen(filename, "w");
	fprintf(fid,"%0.4f\n",Rd);
	fclose(fid);

	printf("%s is done.\n",myname);

    printf("------------------------------------------------------\n");
    now = time(NULL);
    printf("%s\n", ctime(&now));


    free(v);
    free(F);
    free(R);
    free(DetID);
    free(DetW);
    free(DetS);
	free(DetS2);
	free(DetE); // RMT, added probe
    free(DetL);
    free(DetZ);
    return 0;
} /* end of main */



/* SUBROUTINES */

/**************************************************************************
 *	RandomGen
 *      A random number generator that generates uniformly
 *      distributed random numbers between 0 and 1 inclusive.
 *      The algorithm is based on:
 *      W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *      Flannery, "Numerical Recipes in C," Cambridge University
 *      Press, 2nd edition, (1992).
 *      and
 *      D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *      of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *      When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *      When Type is 1, returns a random number.
 *      When Type is 2, gets the status of the generator.
 *      When Type is 3, restores the status of the generator.
 *
 *      The status of the generator is represented by Status[0..56].
 *
 *      Make sure you initialize the seed before you get random
 *      numbers.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

// KE: RandomGen(1,0,NULL): generates a random number between 0 and 1
double RandomGen(char Type, long Seed, long *Status){
    static long i1, i2, ma[56];   /* ma[0] is not used. */
    long        mj, mk;
    short       i, ii;

    if (Type == 0) {              /* set seed. */
        mj = MSEED - (Seed < 0 ? -Seed : Seed);
        mj %= MBIG;
        ma[55] = mj;
        mk = 1;
        for (i = 1; i <= 54; i++) {
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ)
                mk += MBIG;
            mj = ma[ii];
        }
        for (ii = 1; ii <= 4; ii++) // KE: ii is k in random number generator ran3 in mcml
            for (i = 1; i <= 55; i++) {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ)
                    ma[i] += MBIG;
            }
        i1 = 0; // KE: i1 is inext in random number generator ran3 in mcml
        i2 = 31; // KE: i2 is inextp in random number generator ran3 in mcml
    } else if (Type == 1) {       /* get a number. */
        if (++i1 == 56)
            i1 = 1;
        if (++i2 == 56)
            i2 = 1;
        mj = ma[i1] - ma[i2];
        if (mj < MZ)
            mj += MBIG;
        ma[i1] = mj;
        return (mj * FAC);
    } else if (Type == 2) {       /* get status. */
        for (i = 0; i < 55; i++)
            Status[i] = ma[i + 1];
        Status[55] = i1;
        Status[56] = i2;
    } else if (Type == 3) {       /* restore status. */
        for (i = 0; i < 55; i++)
            ma[i + 1] = Status[i];
        i1 = Status[55];
        i2 = Status[56];
    } else
        puts("Wrong parameter to RandomGen().");
    return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


/***********************************************************
 *  Determine if the two position are located in the same voxel
 *	Returns 1 if same voxel, 0 if not same voxel.
 ****/
Boolean SameVoxel(double x1,double y1,double z1, double x2, double y2, double z2,
				double dx,double dy,double dz)
{
    double xmin=min2((floor)(x1/dx),(floor)(x2/dx))*dx;
    double ymin=min2((floor)(y1/dy),(floor)(y2/dy))*dy;
    double zmin=min2((floor)(z1/dz),(floor)(z2/dz))*dz;
    double xmax = xmin+dx;
    double ymax = ymin+dy;
    double zmax = zmin+dz;
    Boolean sv=0;

    sv=(x1<=xmax && x2<=xmax && y1<=ymax && y2<=ymax && z1<zmax && z2<=zmax);
    return (sv);
}

/***********************************************************
 * max2
 ****/
double max2(double a, double b) {
    double m;
    if (a > b)
        m = a;
    else
        m = b;
    return m;
}

/***********************************************************
 * min2
 ****/
double min2(double a, double b) {
    double m;
    if (a >= b)
        m = b;
    else
        m = a;
    return m;
}
/***********************************************************
 * min3
 ****/
double min3(double a, double b, double c) {
    double m;
    if (a <=  min2(b, c))
        m = a;
    else if (b <= min2(a, c))
        m = b;
    else
        m = c;
    return m;
}

/***********************************************************
 *	FRESNEL REFLECTANCE & Refraction - RMT
 * Computes reflectance as photon passes from medium 1 to
 * medium 2 with refractive indices n1,n2. Incident
 * angle a1 is specified by cosine value ca1 = cos(a1).
 * Program returns value of transmitted angle a1 as
 * value in *ca2_Ptr = cos(a2).
 ****/
double ReflectRefraction(double n1, double n2, double* x, double* y, double* z, double* ux, double* uy, double* uz, int* face_dir, int* reflected)
{
	//Determining which axis is the axial one
	double tempux = *ux;
	double tempuy = *uy;
	double tempuz = *uz;
	double ur1;
	double us1;
	double ut1;
	double ur2;
	double us2;
	double ut2;
	double ur22;
	double utot2;
	double cos1;
	double cos2;
	double rs;
	double rp;
	double Rtot;
	
	if(*face_dir == 1) {
		ur1 = *ux;
		us1 = *uy;
		ut1 = *uz;
	} 
	else if(*face_dir == 2) {
		ur1 = *uy;
		us1 = *ux;
		ut1 = *uz;
	} 
	else if(*face_dir == 3) {
		ur1 = *uz;
		us1 = *ux;
		ut1 = *uy;
	}
	
	utot2 = n2/n1; // Not normalize 
	us2 = us1/utot2; // First lateral direction
	ut2 = ut1/utot2; // Second lateral direction
	ur22 = 1-us2*us2-ut2*ut2; //Square of the axial direction
	
	if(ur22 < 0) // In this case, we have a total inter reflection
	{
		ur2 = -ur1; // The light is reflected in the opposite direction
		us2 = us1; // The lateral trajectory doesn't change during reflection
		ut2 = ut1; // The lateral trajectory doesn't change during reflection
		*reflected = 1;
	}
	else
	{
		ur2 = sqrt(ur22); // Axial direction
		if(ur1 < 0) // Make sure to have the proper direction that was lost in the squarre
		{
			ur2 = -ur2; 
		}
		cos1 = ur1; // Cosine of the incident angle
		cos2 = ur2; // Cosine of the transmited angle
		rs = (n1*cos1-n2*cos2)/(n1*cos1+n2*cos2); // Calculating s polarized reflection
		rp = (n1*cos2-n2*cos1)/(n1*cos2+n2*cos1); // Calculating p polarized reflection
		Rtot = 0.5*(rs*rs+rp*rp); // Total reflection, no polarisation implemented yet
		if(RandomNum < Rtot) 
		{
			ur2 = -ur1; // The light is reflected in the opposite direction
			us2 = us1; // The lateral trajectory doesn't change during reflection
			ut2 = ut1; // The lateral trajectory doesn't change during reflection
			*reflected = 1;
		}
		else
		{
			*reflected = 0;
		}
	}
	
	// Calculate the new photon directions
	if(*face_dir == 1) {
		*ux = ur2;
		*uy = us2;
		*uz = ut2;
	} 
	else if(*face_dir == 2) {
		*ux = us2;
		*uy = ur2;
		*uz = ut2;
	} 
	else if(*face_dir == 3) {
		*ux = us2;
		*uy = ut2;
		*uz = ur2;
	}
	
	// Calculate the new photon position
	*x += (-tempux + *ux) * ls; // Take a ls step back to be at the frontiere and take a ls step in the right direction
	*y += (-tempuy + *uy) * ls;
	*z += (-tempuz + *uz) * ls;
	
} /******** END SUBROUTINE **********/


/***********************************************************
 *	Mirror reflection - RMT
 * Gives a reflection with the reflection index chances
 ****/
double MirrorReflection(double n1, double n2, double mr, double* x, double* y, double* z, double* ux, double* uy, double* uz, int* face_dir, int* reflected)
{
	//Determining which axis is the axial one
	double tempux = *ux;
	double tempuy = *uy;
	double tempuz = *uz;
	double ur1;
	double us1;
	double ut1;
	double ur2;
	double us2;
	double ut2;
	double ur22;
	double utot2;
	double cos1;
	double cos2;
	double rs;
	double rp;
	double Rtot;
	
	if(*face_dir == 1) {
		ur1 = *ux;
		us1 = *uy;
		ut1 = *uz;
	} 
	else if(*face_dir == 2) {
		ur1 = *uy;
		us1 = *ux;
		ut1 = *uz;
	} 
	else if(*face_dir == 3) {
		ur1 = *uz;
		us1 = *ux;
		ut1 = *uy;
	}
	
	utot2 = n2/n1; // Not normalize 
	us2 = us1/utot2; // First lateral direction
	ut2 = ut1/utot2; // Second lateral direction
	ur22 = 1-us2*us2-ut2*ut2; //Square of the axial direction
	
	if(ur22 < 0) // In this case, we have a total inter reflection
	{
		ur2 = -ur1; // The light is reflected in the opposite direction
		us2 = us1; // The lateral trajectory doesn't change during reflection
		ut2 = ut1; // The lateral trajectory doesn't change during reflection
		*reflected = 1;
	}
	else
	{
		ur2 = sqrt(ur22); // Axial direction
		if(ur1 < 0) // Make sure to have the proper direction that was lost in the squarre
		{
			ur2 = -ur2; 
		}

		if(RandomNum < mr) 
		{
			ur2 = -ur1; // The light is reflected in the opposite direction
			us2 = us1; // The lateral trajectory doesn't change during reflection
			ut2 = ut1; // The lateral trajectory doesn't change during reflection
			*reflected = 1;
		}
		else
		{
			*reflected = 0;
		}
	}
	
	// Calculate the new photon directions
	if(*face_dir == 1) {
		*ux = ur2;
		*uy = us2;
		*uz = ut2;
	} 
	else if(*face_dir == 2) {
		*ux = us2;
		*uy = ur2;
		*uz = ut2;
	} 
	else if(*face_dir == 3) {
		*ux = us2;
		*uy = ut2;
		*uz = ur2;
	}
	
	// Calculate the new photon position
	*x += (-tempux + *ux) * ls; // Take a ls step back to be at the frontiere and take a ls step in the right direction
	*y += (-tempuy + *uy) * ls;
	*z += (-tempuz + *uz) * ls;
	
} /******** END SUBROUTINE **********/

/***********************************************************
 *	Reverse Erf function - RMT
 * Used to calculate the gaussian beam profil. Is an approximation.
 * Work based on "A handy approximation for the error function and its inverse" by Sergei Winitzki.
 ****/
double myErfInv2(double rndE)
{
	
   double tt1;
   double tt2;
   double lnx;
   double sgn;
   
   sgn = (rndE < 0) ? -1.0f : 1.0f;

   rndE = (1 - rndE)*(1 + rndE);        // x = 1 - x*x;
   lnx = logf(rndE);

   tt1 = 2/(PI*0.147) + 0.5f * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
	
} /******** END SUBROUTINE **********/


/***********************************************************
 *	FRESNEL REFLECTANCE
 * Computes reflectance as photon passes from medium 1 to
 * medium 2 with refractive indices n1,n2. Incident
 * angle a1 is specified by cosine value ca1 = cos(a1).
 * Program returns value of transmitted angle a1 as
 * value in *ca2_Ptr = cos(a2).
 ****/
double RFresnel(double n1,		/* incident refractive index.*/
                double n2,		/* transmit refractive index.*/
                double ca1,		/* cosine of the incident */
                /* angle a1, 0<a1<90 degrees. */
                double *ca2_Ptr) 	/* pointer to the cosine */
/* of the transmission */
/* angle a2, a2>0. */
{
    double r;

    if(n1==n2) { /** matched boundary. **/
        *ca2_Ptr = ca1;
        r = 0.0;
	}
    else if(ca1>(1.0 - 1.0e-12)) { /** normal incidence. **/
        *ca2_Ptr = ca1;
        r = (n2-n1)/(n2+n1);
        r *= r;
	}
    else if(ca1< 1.0e-6)  {	/** very slanted. **/
        *ca2_Ptr = 0.0;
        r = 1.0;
	}
    else  {			  		/** general. **/
        double sa1, sa2; /* sine of incident and transmission angles. */
        double ca2;      /* cosine of transmission angle. */
        sa1 = sqrt(1-ca1*ca1);
        sa2 = n1*sa1/n2;
        if(sa2>=1.0) {
            /* double check for total internal reflection. */
            *ca2_Ptr = 0.0;
            r = 1.0;
		}
        else {
            double cap, cam;	/* cosines of sum ap or diff am of the two */
            /* angles: ap = a1 + a2, am = a1 - a2. */
            double sap, sam;	/* sines. */
            *ca2_Ptr = ca2 = sqrt(1-sa2*sa2);
            cap = ca1*ca2 - sa1*sa2; /* c+ = cc - ss. */
            cam = ca1*ca2 + sa1*sa2; /* c- = cc + ss. */
            sap = sa1*ca2 + ca1*sa2; /* s+ = sc + cs. */
            sam = sa1*ca2 - ca1*sa2; /* s- = sc - cs. */
            r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam);
            /* rearranged for speed. */
		}
	}
    return(r);
} /******** END SUBROUTINE **********/



/***********************************************************
 * the boundary is the face of some voxel
 * find the the photon's hitting position on the nearest face of the voxel
 * and update the step size.
****/
double FindVoxelFace(double x1,double y1,double z1, double x2, double y2, double z2,
					double dx,double dy,double dz, double ux, double uy, double uz)
{
    double x_1 = x1/dx;
    double y_1 = y1/dy;
    double z_1 = z1/dz;
    double x_2 = x2/dx;
    double y_2 = y2/dy;
    double z_2 = z2/dz;
    double fx_1 = floor(x_1) ;
    double fy_1 = floor(y_1) ;
    double fz_1 = floor(z_1) ;
    double fx_2 = floor(x_2) ;
    double fy_2 = floor(y_2) ;
    double fz_2 = floor(z_2) ;
    double x=0, y=0, z=0, x0=0, y0=0, z0=0, s=0;

    if ((fx_1 != fx_2) && (fy_1 != fy_2) && (fz_1 != fz_2) ) { //#10
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        fz_2=fz_1+SIGN(fz_2-fz_1);//added

        x = (max2(fx_1,fx_2)-x_1)/ux;
        y = (max2(fy_1,fy_2)-y_1)/uy;
        z = (max2(fz_1,fz_2)-z_1)/uz;
        if (x == min3(x,y,z)) {
            x0 = max2(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else if (y == min3(x,y,z)) {
            y0 = max2(fy_1,fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
        else {
            z0 = max2(fz_1,fz_2);
            y0 = (z0-z_1)/uz*uy+y_1;
            x0 = (z0-z_1)/uz*ux+x_1;
        }
    }
    else if ( (fx_1 != fx_2) && (fy_1 != fy_2) ) { //#2
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        x = (max2(fx_1,fx_2)-x_1)/ux;
        y = (max2(fy_1,fy_2)-y_1)/uy;
        if (x == min2(x,y)) {
            x0 = max2(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else {
            y0 = max2(fy_1, fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
    }
    else if ( (fy_1 != fy_2) &&(fz_1 != fz_2) ) { //#3
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        y = (max2(fy_1,fy_2)-y_1)/uy;
        z = (max2(fz_1,fz_2)-z_1)/uz;
        if (y == min2(y,z)) {
            y0 = max2(fy_1,fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
        else {
            z0 = max2(fz_1, fz_2);
            x0 = (z0-z_1)/uz*ux+x_1;
            y0 = (z0-z_1)/uz*uy+y_1;
        }
    }
    else if ( (fx_1 != fx_2) && (fz_1 != fz_2) ) { //#4
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        x = (max2(fx_1,fx_2)-x_1)/ux;
        z = (max2(fz_1,fz_2)-z_1)/uz;
        if (x == min2(x,z)) {
            x0 = max2(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else {
            z0 = max2(fz_1, fz_2);
            x0 = (z0-z_1)/uz*ux+x_1;
            y0 = (z0-z_1)/uz*uy+y_1;
        }
    }
    else if (fx_1 != fx_2) { //#5
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        x0 = max2(fx_1,fx_2);
        y0 = (x0-x_1)/ux*uy+y_1;
        z0 = (x0-x_1)/ux*uz+z_1;
    }
    else if (fy_1 != fy_2) { //#6
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        y0 = max2(fy_1, fy_2);
        x0 = (y0-y_1)/uy*ux+x_1;
        z0 = (y0-y_1)/uy*uz+z_1;
    }
    else { //#7
        z0 = max2(fz_1, fz_2);
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        x0 = (z0-z_1)/uz*ux+x_1;
        y0 = (z0-z_1)/uz*uy+y_1;
    }
    //s = ( (x0-fx_1)*dx + (y0-fy_1)*dy + (z0-fz_1)*dz )/3;
    //s = sqrt( SQR((x0-x_1)*dx) + SQR((y0-y_1)*dy) + SQR((z0-z_1)*dz) );
    //s = sqrt(SQR(x0-x_1)*SQR(dx) + SQR(y0-y_1)*SQR(dy) + SQR(z0-z_1)*SQR(dz));
    s = sqrt( SQR((x0-x_1)*dx) + SQR((y0-y_1)*dy) + SQR((z0-z_1)*dz));
    return (s);
}

// KE: the FindVoxelFace2 function is in Listing A4 in Zhao's thesis
// KE: Compute the step size the photon will take to get the first voxel crossing in one single long step.
// KE: We also check whether the photon packet is detected by the assigned detector
/* How much step size will the photon take to get the first voxel crossing in one single
	long step? */ //
double FindVoxelFace2(double x1, double y1, double z1, int* det_num, int Pick_det, double detx, double det_radius, double det_z, double cos_accept, int Ndetectors, double dx, double dy, double dz, double ux, double uy, double uz, int* face_dir)
{

    // KE: ix1, iy1, iz1: indices of the voxel where the photon is currently in
        int ix1 = floor(x1 / dx);
		int iy1 = floor(y1 / dy);
		int iz1 = floor(z1 / dz);
		int izd = floor(det_z / dz);

        // KE: ix2, iy2, iz2: indices of the voxel faces lying ahead of the photon's propagation path
        int ix2, iy2, iz2;
        // KE: equal to equation 4.12 in Zhao's thesis
		if (ux >= 0) ix2 = ix1 + 1;
		else ix2 = ix1;
        // KE: equal to equation 4.13 in Zhao's thesis
		if (uy >= 0) iy2 = iy1 + 1;
		else iy2 = iy1;
        // KE: equal to equation 4.14 in Zhao's thesis
		if (uz >= 0) iz2 = iz1 + 1;
		else iz2 = iz1;
        // KE: xs, ys, zs: distances from these voxel faces to the current position of the photon utilizing its propagation directions
        double xs = fabs((ix2 * dx - x1) / ux); // KE: equal to equations 4.15 in Zhao's thesis
		double ys = fabs((iy2 * dy - y1) / uy); // KE: equal to equations 4.16 in Zhao's thesis
		double zs = fabs((iz2 * dz - z1) / uz); // KE: equal to equations 4.17 in Zhao's thesis
        // KE: s: desired distance of the photon to its closest voxel face
        double s = min3(xs, ys, zs);
		// check detection
		if (-uz >= cos_accept && izd == iz1 && s == zs && fabs(y1 + s * uy) <= det_radius)
        {
			if (fabs(x1 + s * ux - detx) <= det_radius)
                *det_num = Pick_det;
		}
		if (s == xs) *face_dir = 1;		
		else if (s == ys) *face_dir = 2;
		else if (s == zs) *face_dir = 3;
        return (s);
}
double RFresnel(double n1, double n2, double ca1, double *ca2_Ptr);
