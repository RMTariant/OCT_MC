# OCT_MC
Monte-Carlo simulation of an FD-OCT B-Scan.
Add description

how to use

version history
 hello martin hello uday
more stuufff
branch master

Description

make_TissueList_mergedcode_2layers_bloodvessel.m


make_Tissue_mergedcode_2layers_bloodvessel.m

look_mergedcode_2layers_bloodvessel.m

time_min = A(1);
Nx = A(2);
Ny = A(3);
Nz = A(4);
dx = A(5);
dy = A(6);
dz = A(7);
mcflag = A(8);
launchflag = A(9);
boundaryflag = A(10);
xs = A(11);
ys = A(12);
zs = A(13);
xfocus = A(14);
yfocus = A(15);
zfocus = A(16);
ux0 = A(17);
uy0 = A(18);
uz0 = A(19);
radius = A(20);
waist = A(21);
zsurf = A(22);
Nt = A(23);

DetS: Path Length of the detected photons
DetW: Weigth of the detected photons
DetL: Likelihood ratios of detected photons
DetID: IDs of detected photons
F(y,x,z): Data Fluence rate 3D information stored in a 1D vector
count: Number of points of stored information
T(y,x,z): Tissue structure, indexes of the different tissue layers
Tzx: T in the middle of y
R(y,x): Escaping flux at zsurf?
Rd: ???




