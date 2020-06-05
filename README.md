# OCT_MC
Monte-Carlo simulation of an FD-OCT B-Scan.

Description

===Missing Features (to be implemented)===
-Select the number of photons for the simulation
-Fresnel reflection
-Mirror reflection
-Photon pathlength dependent on the refractive index
-Refraction
-Noise
-Gaussian beam source
-Dispersion?
-Wavelength width on the detection
-Focusing/diverging beams instead of perfectly foccus
-Check cone of acceptance in the fiber

MergedCode.c

===Functions===


===Variables===
a: Variable used in the original code to determine if two positions are in the same voxel
b: Variable used in the original code to determine if two positions are in the same voxel
c: Variable used in the original code to determine if two positions are in the same voxel
W: Photon weigth (max(W) = 1)
photon_status: ALIVE or DEAD?
i_photon: increment of photon photon packet
c_photon: count of collected photon packet
j: increment of a for loop. used twice?
NN: Number of voxel in the 3D simulated volume
Nyx: Number of voxel in a 2D cut on the y and x axis of the simulated volume
Nx,Ny,Nz: Size of the simualted space in voxel in the corresponding direction
dx,dy,dz: Width of a voxel in cm
xs,ys,zs: Control of the position of the source in cm? xs and ys = 0 for a centered source.
x,y,z: Current position of the photon
x_cont,y_cont,z_cont: x,y,z coordonates of the continuing photon that "died"?
tempx,tempy,tempz: Temporary value for x,y,z.
ux,uy,uz: Direction of the photon
ux_cont,uy_cont,uz_cont: ux,uy,uz coordonates of the continuing photon that "died"?
ix,iy,iz: Indexes used to determine if the photon goes through a border and change medium
detx: Position of the dectector/launch position. Ranges from -radius to +radius
det_z: Position of the dectector in z. Redundant, always = zs.
F: Fluence rate
R: Reflectance It was specify in the original code that it was not ready to be used. We don't use it either.
Rd: I think it quantifies the total quantity of photon reflected at the surface. Not used.
det_num: Detector ID. Can have a value of -1 if no detection has occured yet.
first_bias_done: 0 or 1? Determine if the photon has been bias or not yet.
cont_exist: 0 or 1? determine if the photon have been split yet.
L_current: Photon likelihood. Initially = 1
L_cont: L_current of the continuing photon that "died"?
s_total: photon total path length. Starts at 0.
s_total_cont: s_total of the continuing photon that "died"?
z_max: photon depth reached
z_max_cont: z_max of the continuing photon that "died"?
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
boundaryflag: Determines where the photon can escape. 0 = no boundaries, 1 = escape at all boundaries, 2 = escape at surface only From the original code.
rnd: Random variable ranging from 0 to 1 sometimes including or excluding one or both of the extremities
sleft: Distance the photon have left to propagate but without considering mus. (Better formulate)
s: Random propagating distance in cm
sv: 1 = The photon is in the same voxel as before. 0 = The photon have changed voxel.
absorb: Proportion of the photon packet absorb
W: Weigth of the photon
W_cont: Weigth of the photon of the continuing photon that "died"?
i: Index of the voxel. All voxels in the simulated space have a specific index associated to it.
i_cont: i of the continuing photon that "died"?
ls: A tiny step value to ensure that the photon is not exactly between two voxels. In cm.
v: Related to the tissue type
type: Current tissue type where the photon is located
mua: Current absorption coeficient in cm
muav: Possible absorption coeficients in the simulation in cm
mus: Current scattering coeficient in cm
musv: Possible scattering coeficients in the simulation in cm
g: Current anisotropy
gv: Possible anisotropy in the simulation
cont_exist: Determines if a photon split will occures 0=no and 1=yes
first_bias_done: Determines if a scaterred bias has occured yet =1 or not =0


DetS: Path Length of the detected photons
DetW: Weigth of the detected photons
DetL: Likelihood ratios of detected photons
DetID: IDs of detected photons (range from 1 to 512 for 512 detectors/emittors)
c_photon: Current photon packet itteration?




