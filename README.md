# OCT_MC
Monte-Carlo simulation of an FD-OCT B-Scan.

Description

Variables


F(y,x,z): Data Fluence rate 3D information stored in a 1D vector
count: Number of points of stored information
T(y,x,z): Tissue structure, indexes of the different tissue layers
Tzx: T in the middle of y
R(y,x): Escaping flux at zsurf?
Rd: ???

MergedCode.c

a: Variable used in the original code to determine if two positions are in the same voxel
b: Variable used in the original code to determine if two positions are in the same voxel
c: Variable used in the original code to determine if two positions are in the same voxel
W: Photon weigth (max(W) = 1)
photon_status: ALIVE or DEAD?
i_photon: increment of photon photon packet
c_photon: count of collected photon packet (by detector I guess)
j: increment of a for loop. used twice?
NN: Number of voxel in the 3D simulated volume
Nyx: Number of voxel in a 2D cut on the y and x axis of the simulated volume
Nx,Ny,Nz: Size of the simualted space in voxel in the corresponding direction
dx,dy,dz: Width of a voxel in cm
xs,ys,zs: Control of the position of the source in cm? xs and ys = 0 for a centered source.
x,y,z: Current position of the photon
tempx,tempy,tempz: Temporary value for x,y,z.
ux,uy,uz: Initial direction of the photon
ix,iy,iz: Indexes used to determine if the photon goes through a border and change medium
detx: Position of the dectector/launch position. Ranges from -radius to +radius
det_z: Position of the dectector in z. Redundant, always = zs.
F: Fluence rate
R: Reflectance It was specify in the original code that it was not ready to be used. We don't use it either.
det_num: Detector ID. Can have a value of -1 if no detection has occured yet.
first_bias_done: 0 or 1? Determine if the photon has been bias or not yet.
cont_exist: 0 or 1? determine if the photon have been split yet.
L_current: Photon likelihood. Initially = 1
s_total: photon total path length. Starts at 0.
z_max: photon depth reached
Ndetectors: Number of detectors/Alines
det_radius: Radius of the detector in cm
cos_accept: Aperture of the detector/cos of the angle
Pick_det: Picked detector. Ranged from 1 to Ndetectors in 
radius: NOT A RADIUS. It is half the B-line width.
detx: Position in cm? of the current A-line position
launchflag: if = 1 Manually set launch direction to straigth down.
mcflag: 0 = uniform beam. Other values are for 1 = gaussian, 2 = isotropic and 3 rectangular beam
bflag: 1 = The photon is still inside the simulated space. 0 = it no longer is.
surfflag: 0 = photon is inside the tissue part of the simulated space. 1 = photon is inside the air part of the simulated space.
rnd: Random variable ranging from 0 to 1 sometimes including or excluding one or both of the extremities
sleft: Distance the photon have left to propagate but without considering mus. (Better formulate)
s: Random propagating distance in cm
sv: 1 = The photon is in the same voxel as before. 0 = The photon have changed voxel.
absorb: Proportion of the photon packet absorb
W: Weigth of the photon
i: Index of the voxel. All voxels in the simulated space have a specific index associated to it.
ls: A tiny step value to ensure that the photon is not exactly between two voxels. In cm.

DetS: Path Length of the detected photons
DetW: Weigth of the detected photons
DetL: Likelihood ratios of detected photons
DetID: IDs of detected photons (range from 1 to 512 for 512 detectors/emittors)
c_photon: Current photon packet itteration?
s_total: Total path length of the current photon packet
s: Path length of the current scattering event
surfflag: 0: photon inside tissue, 1: photon outside tissue?
bflag: Boundarty flag. Determines if the photon is still inside the (simulated?) volume (1) or not (0)
boundaryflag: 0 = no boundary, 1 = escape at all boundaries, 2 = escape at surface only




