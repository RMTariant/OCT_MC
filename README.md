# OCT_MC
Monte-Carlo simulation of an FD-OCT B-Scan.

Description

Variable


F(y,x,z): Data Fluence rate 3D information stored in a 1D vector
count: Number of points of stored information
T(y,x,z): Tissue structure, indexes of the different tissue layers
Tzx: T in the middle of y
R(y,x): Escaping flux at zsurf?
Rd: ???

MergedCode.c
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
mcflag: 0 = uniform flat field, 1 = Gaussian beam, 2 = isotropic point source, 3 = squarre source
launchflag: 0 = ignore, 1 = manually set
boundaryflag
ux0: Initial launch direction?
uy0: Initial launch direction?
uz0: Initial launch direction?
x, y, z: position?
xs, ys, zs: posution?
detx: Detector x position?





