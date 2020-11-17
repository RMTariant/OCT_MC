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
F = F/sum(F(:));
figure
x = (0:Nx-1)*dx;
z = (0:Nz-1)*dz;
imagesc(x,z,squeeze(db(F(round(Ny/2),:,:))));
c = colorbar;
c.Label.String = 'decibel';

%% other stuff








