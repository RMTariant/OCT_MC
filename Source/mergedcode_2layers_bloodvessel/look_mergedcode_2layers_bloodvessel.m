home; clear
format compact
commandwindow

SAVEPICSON = 1;
if SAVEPICSON
    sz = 10; fz = 7; fz2 = 5; % to use savepic.m
else
    sz = 12; fz = 9; fz2 = 7; % for screen display
end

myname = 'skinvessel2layersb'; nm = 532;

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%% parameters
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
j = 23;
for i=1:Nt
    j=j+1;
    muav(i,1) = A(j);
    j=j+1;
    musv(i,1) = A(j);
    j=j+1;
    gv(i,1) = A(j);
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
    DetID = fread(fid, 'float');
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
x = ([1:Nx]-Nx/2-1/2)*dx; %([1:200] - 100 - 0.5) * 4*10^-4; [-0.0398,0.0398] with a difference of 4*10^-4
y = ([1:Ny]-Ny/2-1/2)*dy; %([1:200] - 100 - 0.5) * 4*10^-4; [-0.0398,0.0398] with a difference of 4*10^-4
z = ([1:Nz]-1/2)*dz; %([1:200]-1/2)*4*10^-4; [2*10^-4, 0.0798] with a difference of 4*10^-4
ux = [2:Nx-1]; %[2:199]
uy = [2:Ny-1]; %[2:199]
uz = [2:Nz-1]; %[2:199]
zmin = min(z); % 2*10^-4
zmax = max(z); % 0.0798
zdiff = zmax-zmin; % 0.0796
xmin = min(x); % -0.0398
xmax = max(x); % 0.0398
xdiff = xmax-xmin; % 0.0796

%% Look at structure, Tzx
Tzx = reshape(T(Ny/2,:,:),Nx,Nz)'; % reshape(T(100,200,200)) 200x200 double
tissue = make_TissueList_mergedcode_2layers_bloodvessel(nm);
Nt = length(tissue); % different types of tissue, in this case:2

figure(1);clf
imagesc(x(ux),z(uz),Tzx(uz,ux),[1 Nt])
hold on
cmap = makecmap(Nt);
colormap(cmap)
colorbar
set(gca,'fontsize',sz)
set(colorbar,'fontsize',1)
xlabel('x [cm]')
ylabel('z [cm]')
title('Tissue','fontweight','normal','fontsize',fz2)
for i=1:Nt
    yy = zmin + (Nt-i)/(Nt-1)*zdiff;
    text(xmax*1.4,yy, sprintf('%d %s',i,tissue(i).name),'fontsize',fz2)
end

% draw launch
N = 10; % # of beam rays drawn
switch mcflag
    case 0 % uniform
        for i=0:N
            plot((-radius + 2*radius*i/N)*[1 1],[zs max(z)],'r-')
        end

    case 1 % Gaussian
        for i=0:N
            plot([(-radius + 2*radius*i/N) xfocus],[zs zfocus],'r-')
        end

    case 2 % iso-point
        for i=1:N
            th = (i-1)/19*2*pi;
            xx = Nx/2*cos(th) + xs;
            zz = Nx/2*sin(th) + zs;
            plot([xs xx],[zs zz],'r-')
        end
        
    case 3 % rectangle
        zz = max(z);
        for i=0:N
            xx = -radius + 2*radius*i/N;
            plot([xx xx],[zs zz],'r-')
        end
end
axis equal image

if SAVEPICSON
    name = sprintf('%s_tissue.fig',myname);
    %savefig(1,[4 3],name)
    savefig(name)
end


%% Look at Fluence Fzx @ launch point
Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source
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
colormap(makec2f)
%colormap(map) sets the colormap for the current figure to the colormap
%specified by map
axis equal image
%axis([min(x) max(x) min(z) max(z)])
text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',time_min),...
    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%s_Fzx.fig',myname);
    savefig(name)
end

%% look Fzy
Fzy = reshape(F(:,Nx/2,:),Ny,Nz)';

iy = round((dy*Ny/2 + 0.15)/dy);
iz = round(zs/dz);
zzs  = zs;
%Fdet = mean(reshape(Fzy(iz+[-1:1],iy+[0 1]),6,1));

figure(3);clf
imagesc(y,z,log10(Fzy))
hold on
text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('y [cm]')
ylabel('z [cm]')
title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
colormap(makec2f)
axis equal image
text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',time_min),...
    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%s_Fzy',myname);
    savefig(name)
end

%% look Azx
Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source
mua = muav(reshape(T(Ny/2,:,:),Nx,Nz)');
Azx = Fzx.*mua;

figure(4);clf
imagesc(x,z,log10(Azx))
hold on
text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( A )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('Deposition A [W/cm^3/W.delivered] ','fontweight','normal','fontsize',fz)
colormap(makec2f)
axis equal image
%axis([min(x) max(x) min(z) max(z)])
text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',time_min),...
    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%s_Azx',myname);
    savefig(name)
end

%%
%{
figure(5);clf
ux = [2:Nx-1]';
U = Ryx(ux,ux);
Rtot = sum(U(:))*dx*dy

imagesc(x(ux),y(ux),U,[0 mean(U(:))*2])
colorbar
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('y [cm]')
s = sprintf('R [W/cm^2/W.delivered]\nsum(R(y,x))dxdy = %0.3f, R_d = %0.3f',...
    Rtot,Rd);
title(s)
axis equal image

figure(6);clf
plot(x(ux),U(Ny/2,:),'r.-')
hold on


drawnow

cd(homedir)
ifigs=1:5;
disp('done')
%}


