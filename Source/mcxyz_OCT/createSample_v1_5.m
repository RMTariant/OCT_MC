% Create the files for the simulation
   
clear
format compact
clc
home


%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci 
                        % 0 = don't save. Just check the program.

myname      = 'Trial1';% name for files: myname_T.bin, myname_H.mci  
time_min    = 0.01;      	% time duration of the simulation [min] <----- run time -----
Nx          = 100;    	% # of bins in each dimension of cube 
Ny          = 100;    	% # of bins in each dimension of cube 
Nz          = 100;    	% # of bins in each dimension of cube 
binsize     = 0.003;     	% size of each bin, eg. [cm]

% Set Monte Carlo launch flags (not in use)
mcflag      = 0;     	% launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 
                        % 3 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
launchflag  = 1;        % 0 = let mcxyz.c calculate launch trajectory
                        % 1 = manually set launch vector.
boundaryflag = 1;       % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x, y, bottom z
                        % boundaries

% Sets position of source with 0 centered on each Aline
xs          = 0;      	% x of source
ys          = 0;        % y of source
zs          = 1e-8; % z of source must start in simulation

% Set position of focus, so mcxyz can calculate launch trajectory (not in use)
xfocus      = 0;        % set x,position of focus
yfocus      = 0;        % set y,position of focus
zfocus      = inf;    	% set z,position of focus (=inf for collimated beam)

% Set detection parameter
radius      = 0.05;     % Half width of the BScan
waist       = 0.05;  	% Width of the scanned beam (Not in use)
Ndetectors  = 800;      % Number of Aline per BScan
flens       = 0.06;        %Focal lens in m 0.01 or 0.06
flens = 40e-3;
beamw       = 0.002;        %Beam diameter at imaging lens in m
beamw = 35e-5;
det_radius  = 1310e-7*2/pi/atan(beamw/2/flens);    % Width of the beam at the imaging lens
cos_accept  = cos(atan(beamw/2/flens)); % Cos of the accepted angle

% Parameters for gaussian beam 
lambda = 1310e-9;
f = flens;
D = 35e-5;
z_f_img = -20e-3;
h_step = 1e-4; 

% only used if launchflag == 1 (manually set launch trajectory): (not in use)
ux0         = 0;      % trajectory projected onto x axis
uy0         = 0;      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1

% Bias scattering parameter
a_coef    = 0.9;
p         = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
% Create Sample
%%%
zsurf = 0.0;  % position of air/skin surface
T = double(zeros(Nx,Ny,Nz));
T(:) = 1;
Nt = 1; %Number if layers
muav(1) = 1.472e-1;    %Absorption coef. of first layer
% musv(1) = 10;   %Scattering coef. of first layer
gv(1) = 0.5;     %Anisotropy of first layer
nrv(1) = 2;      %Refraction index of first layer
%musv(1) = 6.58e0/(1-gv(1));       %Change reduced scattering in scattering
musv(1) = 1;
% muav(2) = 0.1;    %Absorption coef. of first layer
% musv(2) = 1000;   %Scattering coef. of first layer
% gv(2) = 0.7;     %Anisotropy of first layer
% nrv(2) = 1;      %Refraction index of first layer
% T(:,:,2) = 2;

%%%%%%%%%% 
% Prepare Monte Carlo 
%%%


% Specify Monte Carlo parameters    
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



%%
if SAVEON
    tic
    % convert T to linear array of integer values, v(i)i = 0;
    v = uint8(reshape(T,Ny*Nx*Nz,1));

    %% WRITE FILES
    % Write myname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    disp(sprintf('--------create %s --------',myname))
    filename = sprintf('%s_H.mci',myname);
    fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.2f\n',time_min);
        fprintf(fid,'%0.4f\n',a_coef);
        fprintf(fid,'%0.4f\n',p);
        fprintf(fid,'%0.4f\n',Ndetectors);
        fprintf(fid,'%0.6f\n',det_radius);
        fprintf(fid,'%0.6f\n',cos_accept);
        fprintf(fid,'%0.4f\n',lambda);
        fprintf(fid,'%0.4f\n',f);
        fprintf(fid,'%0.4f\n',D);
        fprintf(fid,'%0.4f\n',z_f_img);
        fprintf(fid,'%0.4f\n',h_step);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.4f\n',dx);
        fprintf(fid,'%0.4f\n',dy);
        fprintf(fid,'%0.4f\n',dz);
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
        fprintf(fid,'%0.4f\n',zsurf);
        % tissue optical properties
        fprintf(fid,'%d\n',Nt);
        for i=1:Nt
            fprintf(fid,'%0.6f\n',muav(i));
            fprintf(fid,'%0.6f\n',musv(i));
            fprintf(fid,'%0.6f\n',gv(i));
            fprintf(fid,'%0.6f\n',nrv(i));
        end
    fclose(fid);

    %% write myname_T.bin file
    filename = sprintf('%s_T.bin',myname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,v,'uint8');
    fclose(fid);

    toc
end % SAVEON


%% Look at structure of Tzx at iy=Ny/2
Txzy = shiftdim(T,1);   % Tyxz --> Txzy
Tzx  = Txzy(:,:,round(Ny/2))'; % Tzx

%%
figure; clf
Tzx = squeeze(T(:,round(Ny/2),:))';
imagesc(x,z,Tzx,[0 Nt])
xlabel('x [cm]')
ylabel('z [cm]')
title('\rm Plane y=0')
colorbar
