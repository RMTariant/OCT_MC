% home; clear
% format compact
% commandwindow
% 
% SAVEPICSON = 1;
% if SAVEPICSON
%     sz = 10; fz = 7; fz2 = 5; % to use savepic.m
% else
%     sz = 12; fz = 9; fz2 = 7; % for screen display
% end

myname = 'fluencecompared'; nm = 532;
cd('C:\Users\raphi\Documents\Doctorat\Uday simulation\data')

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%% parameters
Nphoton = A(1);
Nx = A(5);
Ny = A(6);
Nz = A(7);
dx = A(8);
dy = A(9);
dz = A(10);
mcflag = A(11);
launchflag = A(12);
boundaryflag = A(13);
xs = A(14);
ys = A(15);
zs = A(16);
xfocus = A(17);
yfocus = A(18);
zfocus = A(19);
ux0 = A(20);
uy0 = A(21);
uz0 = A(22);
radius = A(23);
waist = A(24);
zsurf = A(25);
Nt = A(26);
j = 26;
for i=1:Nt
    j=j+1;
    muav(i,1) = A(j);
    j=j+1;
    musv(i,1) = A(j);
    j=j+1;
    gv(i,1) = A(j);
    j=j+1;
    nr(i,1) = A(j);
    j=j+1;
end

reportHmci(myname)

%% Load path lengths of detected photons DetS
filename = sprintf('%s_DetS.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetS = fread(fid, 'float');
    fclose(fid);
toc

%% Load path lengths of detected photons DetS
filename = sprintf('%s_DetS2.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetS2 = fread(fid, 'float');
    fclose(fid);
toc

%% Load weights of detected photons DetW
filename = sprintf('%s_DetW.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetW = fread(fid, 'float');
    fclose(fid);
toc

%% Load likelihood ratios of detected photons DetL
filename = sprintf('%s_DetL.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetL = fread(fid, 'float');
    fclose(fid);
toc

%% Load IDs of detected photons DetID
filename = sprintf('%s_DetID.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetID = fread(fid, 'int');
    fclose(fid);
toc

%% Load E probing varaible
filename = sprintf('%s_DetE.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetE = fread(fid, 'float');
    fclose(fid);
toc

%% Load Fluence rate F(y,x,z) 
filename = sprintf('%s_F.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)

%%
% Load tissue structure in voxels, T(y,x,z) 
filename = sprintf('%s_T.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'uint8');
    fclose(fid);
toc
T = reshape(Data,Ny,Nx,Nz); % T(y,x,z)

clear Data

%% Load Escaping flux R(y,x) at zsurf
filename = sprintf('%s_Ryx.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx, 'float');
    fclose(fid);
toc
Ryx = reshape(Data,Ny,Nx); % Ryx(y,x)

%% Rd
filename = sprintf('%s_Rd.dat',myname);
disp(['loading ' filename])
Rd = load(filename);
fprintf('Rd = %0.4f\n',Rd)

%%
F = F/sum(F(~isnan(F(:))));
figure
x = (0:Nx-1)*dx;
z = (0:Nz-1)*dz;
Fzx = reshape(F(round(Ny/2),:,:),Nx,Nz)'; % in z,x plane through source
imagesc(x,z,log10(Fzx));
c = colorbar;
c.Label.String = 'decibel';

%% Look at Fluence Fzx @ launch point
Fzx = reshape(F(round(Ny/2),:,:),Nx,Nz)'; % in z,x plane through source
%Fzx_norm = Fzx - min(Fzx(:));
%Fzx_norm = Fzx ./ max(Fzx(:));
figure(2);clf %deletes from the current figure all graphic objects
imagesc(x,z,log10(Fzx)) % specifies image location
% x and z specify the locations of the corner corresponding to 
% log10(Fzx)(1,1) and log10(Fzx)(m,n)
%imagesc(x,z,log10(Fzx_norm))
% displays the data in array log10(Fzx)as an image that uses the full range
% of colors in the colormap. Each element of log10(Fzx) specifies the color of one pixel of the
% image. The resulting image is an m-by-n grid of pixels with m rows and n
% columns in log10(Fzx). The row and column indices of the elements
% determine the centers of the corresponding pixels.
hold on
text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
colorbar % displays a vertical colorbar to thhe right of the current axes or chart
% displays the current colormap and indicate the mapping of data values
% into the colormap
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
%colormap(makec2f)
%colormap(map) %sets the colormap for the current figure to the colormap
%specified by map
axis equal image
%axis([min(x) max(x) min(z) max(z)])
%text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',time_min),...
%    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%s_Fzx.fig',myname);
    savefig(name)
end

%% Look at attenuation at center

att2 = squeeze(fluence.data(:,51,51));
att1 = squeeze(F(51,51,:));

att22 = att2(40:70);
att12 = att1(40:70);

