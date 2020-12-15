% RMT
% Description: This creates the tissue and simualtor parameters
%% Simulation Parameters


clear
format compact
clc
home

%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci 
                        % 0 = don't save. Just check the program.

% Set Monte-Carlo environment size and time of simulation
cd('C:\Users\raphi\Documents\Doctorat\Uday simulation\data') % RMT files location
myname      = 'tissue';% Name of the saved file
time_min    = 1;      	% Time of the simulation if the number of photon is not specified.
Nphotons    = 1e7;      % Number of photons used in the simulation. Overwrites the time of simulation
Nbins       = 2;    	% # of bins in each dimension of cube Nbins     
binsize     = 1; 	    % size of each bin, eg. [mm] 

% Bias scattering parameters
a_coef      = 0.9;      % Biasing coefficient of the importance sampling
p           = 0.5;      % Probability of additional bias scattering.

% Detector parameters. The detector correspond to the beam at the exit of
% the imaging lens.
N_det       = 800;      % Number of detector in the simulation. They are positioned on a line Span_det long equaly spaced. It corresponds to the number of Alines in a Bscan.
span_det    = 1;        % Number of 
cos_accept  = 0.99619;  % Cosine of the accepted angle of photon being detected Bugged. Need to be added in the .c
det_radius  = 0.001;    % Radius of the detector [cm] at the detector (warning, in the case of OCT, this would be the beam radius at the imaging lens) Need to be added in the .c


% Set Monte Carlo launch flags
mcflag      = 1;     	% launch: 0 = uniform beam, 1 = Gaussian
launchflag  = 0;        % 0 = let mcxyz.c calculate launch trajectory
                        % 1 = manually set launch vector.
boundaryflag = 2;       % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x, y, bottom z
                        % boundaries

% Sets position of source
xs          = 0;      	% x of source
ys          = 0;        % y of source
zs          = 0.0001;  	% z of source zs  Must be inside the simulation

% Set position of focus, so mcxyz can calculate launch trajectory
xfocus      = 0;        % set x,position of focus
yfocus      = 0;        % set y,position of focus
zfocus      = inf;    	% set z,position of focus (=inf for collimated beam)

% only used if mcflag == 0 or 1 or 3 (not 2=isotropic pt.)
radius      = 0.0100;   % 1/e radius of beam at tissue surface
waist       = 2.08e-3;  % 1/e2 radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory):
ux0         = 0;      % trajectory projected onto x axis
uy0         = 0;      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1

%% Prepare Monte Carlo 


Nx = Nbins;
Ny = Nbins;
Nz = Nbins;
dx = binsize;
dy = binsize;
dz = binsize;
x  = ([1:Nx]'-Nx/2)*dx;
y  = ([1:Ny]'-Ny/2)*dy;
z  = [1:Nz]'*dz;
zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);

if isinf(zfocus), zfocus = 1e12; end

zsurf = 0.000; %position of air/skin surface

%% Create tissue

Nt = 3;

%Optical properties, absorption, dispersion and anisotropy. [mm]
muav(1)  = 0.0001;
musv(1)  = 0.0001;
gv(1)    = 0.99;
muav(2)  = 0.0577;
musv(2)  = 0.2965;
gv(2)    = 0.5;
muav(3)  = 1;
musv(3)  = 200;
gv(3)   = 0.8;

% Adding reflection chances
mr(1) = 0;
mr(2) = 0;
mr(3) = 0;

% Creating the tissue shape
T = double ( zeros (Ny ,Nx ,Nz));

%Assigning the tissue type
T(:) = 2;
T(:,:,1) = 1;

 
%% WRITE FILES
if SAVEON
    
    % convert T to linear array of integer values, v(i)i = 0;
    v = uint8(reshape(T,Ny*Nx*Nz,1));
    % Write myname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    disp(sprintf('--------create %s --------',myname))
    filename = sprintf('%s_H.mci',myname);
    fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.0f\n',time_min);
        fprintf(fid,'%0.0f\n',Nphotons);
        fprintf(fid,'%0.4f\n',a_coef);
        fprintf(fid,'%0.4f\n',p);
        fprintf(fid,'%0.0f\n',N_det);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.6f\n',dx);
        fprintf(fid,'%0.6f\n',dy);
        fprintf(fid,'%0.6f\n',dz);
        % launch parameters
        fprintf(fid,'%d\n'   ,mcflag);
        fprintf(fid,'%d\n'   ,launchflag);
        fprintf(fid,'%d\n'   ,boundaryflag);
        fprintf(fid,'%0.4f\n',xs);
        fprintf(fid,'%0.4f\n',ys);
        fprintf(fid,'%0.4f\n',zs);
        fprintf(fid,'%0.4f\n',xfocus);
        fprintf(fid,'%0.4f\n',yfocus);
        fprintf(fid,'%0.4f\n',zfocus);
        fprintf(fid,'%0.4f\n',ux0); % if manually setting ux,uy,uz
        fprintf(fid,'%0.4f\n',uy0);
        fprintf(fid,'%0.4f\n',uz0);
        fprintf(fid,'%0.4f\n',radius);
        fprintf(fid,'%0.4f\n',waist);
        fprintf(fid,'%0.4f\n',zsurf)
        % tissue optical properties
        fprintf(fid,'%d\n',Nt);
        for i=1:Nt
            fprintf(fid,'%0.4f\n',muav(i));
            fprintf(fid,'%0.4f\n',musv(i));
            fprintf(fid,'%0.4f\n',gv(i));
            %fprintf(fid,'%0.4f\n',nr(i));
            %fprintf(fid,'%0.4f\n',mr(i));
        end
    fclose(fid);

    % write myname_T.bin file
    filename = sprintf('%s_T.bin',myname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,v,'uint8');
    fclose(fid);

    
end % SAVEON


%% Look at structure of Tzx at iy=Ny/2
Txzy = shiftdim(T,1);   % Tyxz --> Txzy
Tzx  = Txzy(:,:,round(Ny/2))'; % Tzx

figure(2); clf
sz = 12;  fz = 10; 
imagesc(x,z,Tzx,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('\rm # of A-Scan at y=0: 512')
colorbar
cmap = makecmap(Nt);
colormap(cmap)
set(colorbar,'fontsize',1)
% label colorbar
zdiff = zmax-zmin;

disp('done')

