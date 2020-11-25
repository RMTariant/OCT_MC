% make_Tissue_mergedcode_2layers_bloodvessel.m
% myname = skinvessel_mergedcode_2layers_b % name for files

%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       myname_T.bin = a tissue structure file
%  
%   Also prepares a listing of the optical properties at chosen wavelength
%   [mua, mus, g], for each tissue type specified in myname_T.bin. 
%   This listing is saved in
%       myname_H.mci = the input file for use c code
%
%   Will generate a figure illustrating the tissue with its various
%   tissue types and the beam being launched.
%
%   Uses
%       make_TissueList_mergedcode_2layers_bloodvessel.m
%
%   To use, 
%       1. Prepare make_TissueList_mergedcode_2layers_bloodvessel.m
%       so that it contains the tissue types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, make_Tissue_mergedcode_2layers_bloodvessel.m
%
%   

clear
format compact
clc
home

%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci 
                        % 0 = don't save. Just check the program.

cd('C:\Users\raphi\Documents\Doctorat\Uday simulation\data') % RMT files location
myname      = 'fluencea09p02';%_2layers_bloodvessel';% name for files: myname_T.bin, myname_H.mci  
time_min    = 1;      	% RMT No longuer used time duration of the simulation [min] <----- run time -----
Nphotons    = 2e6;      % RMT Number of photons used in the simulation.
a_coef      = 0.9;      % RMT Biasing coefficient of the importance sampling
p           = 0.2;      % RMT Probability of additional bias scattering.
Ndetectors  = 1;        % RMT Number of detector in the simulation
nm          = 532;   	% desired wavelength of simulation
Nbins       = 101;    	% # of bins in each dimension of cube  Nbins     
binsize     = 0.0001; 	% size of each bin, eg. [cm] or [mm] binsize   

% Set Monte Carlo launch flags
mcflag      = 1;     	% launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 
                        % 3 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
launchflag  = 1;        % 0 = let mcxyz.c calculate launch trajectory
                        % 1 = manually set launch vector.
boundaryflag = 2;       % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x, y, bottom z
                        % boundaries

% Sets position of source
xs          = 0;      	% x of source
ys          = 0;        % y of source
zs          = 0.0001;  	    % z of source zs  

% Set position of focus, so mcxyz can calculate launch trajectory
xfocus      = 0;        % set x,position of focus
yfocus      = 0;        % set y,position of focus
zfocus      = inf;    	% set z,position of focus (=inf for collimated beam)

% only used if mcflag == 0 or 1 or 3 (not 2=isotropic pt.)
radius      = 0.0101;   % 1/e radius of beam at tissue surface
waist       = 0.002;  	% 1/e2 radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory):
ux0         = 0.7;      % trajectory projected onto x axis
uy0         = 0.4;      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% 
% Prepare Monte Carlo 
%%%

Nt = 3;

muav(1)  = 0.001;
musv(1)  = 0.01;
gv(1)    = 0.9;
muav(2)  = 5;
musv(2)  = 1000;
gv(2)    = 0.9;
muav(3)  = 10;
musv(3)  = 2000;
gv(3)   = 0.8;

% Adding refractive index
nr(1) = 1.3; %%%Edit the value
nr(2) = 1.3; %%%Edit the value
nr(3) = 1.4; %%%Edit the value

% Adding reflection chances
mr(1) = 0;
mr(2) = 0;
mr(3) = 0;

% Specify Monte Carlo parameters    
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

%%%%%%
% CREATE TISSUE STRUCTURE T(y,x,z)
%   Create T(y,x,z) by specifying a tissue type (an integer)
%   for each voxel in T.
%
%   Note: one need not use every tissue type in the tissue list.
%   The tissue list is a library of possible tissue types.

T = double ( zeros (Ny ,Nx ,Nz));
T(:) = 2;
T(:,:,1:2) = 1;
T(:,:,31:end) = 3;

zsurf = 0.000; %position of air/skin surface

% %Create three layers
% T(:,:,1:40) = 1;
% T(:,:,41:100) = 2;
% T(:,:,101:200) = 3;
% clear T
% T(1:106,1:106,1:20) = 1;
% T(1:106,1:106,21:66) = 2;
% T(1:106,1:106,66:106) = 3;
% T(1:106,1:106,66) = 4;
 
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
        fprintf(fid,'%0.0f\n',Nphotons);
        fprintf(fid,'%0.4f\n',a_coef);
        fprintf(fid,'%0.4f\n',p);
        fprintf(fid,'%0.0f\n',Ndetectors);
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
            fprintf(fid,'%0.4f\n',nr(i));
            fprintf(fid,'%0.4f\n',mr(i));
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
%%%
% 
% for i=1:Nt
%     yy = (Nt-i)/(Nt-1)*Nz*dz;
%     text(max(x)*1.2,yy, tissue(i).name,'fontsize',fz)
% end
% 
% text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
% axis equal image
% axis([xmin xmax zmin zmax])
% 
% %%% draw launch
% N = 10; % # of beam rays drawn
% switch mcflag
%     case 0 % uniform
%         for i=0:N
%             plot((-radius + 2*radius*i/N)*[1 1],[zs max(z)],'r-')
%         end
% 
%     case 1 % Gaussian
%         for i=0:N
%             plot([(-radius + 2*radius*i/N) xfocus],[zs zfocus],'r-')
%         end
% 
%     case 2 % iso-point
%         for i=1:N
%             th = (i-1)/19*2*pi;
%             xx = Nx/2*cos(th) + xs;
%             zz = Nx/2*sin(th) + zs;
%             plot([xs xx],[zs zz],'r-')
%         end
%         
%     case 3 % rectangle
%         zz = max(z);
%         for i=1:N
%             xx = -radius + 2*radius*i/20;
%             plot([xx xx],[zs zz],'r-')
%         end
% end

%savefig(strcat(myname, '.fig'));
disp('done')

