% lookmcxyz.m
%   Looks at myname_F.bin, created by mcxyz.c 
%   where myname is the name of the run: myname_T.bin, myname_H.mci
%   Makes figures:
%       myname_tissue.jpg   = tissue structure (shows tissue types)
%       myname_Fzx.jpg      = fluence rate vs z,x
%       myname_Fzy.jpg      = fluence rate vs z,y
%   Uses:
%       myname_H.mci    = input file from maketissue.m
%       myname_T.bin    = tissue input file from maketissue.m
%       myname_F.bin    = fluence rate output from Monte Carlo
%       reportH_mci.m   = lists input parameters in myname_H.mci
%       makecmap.m      = makes colormap for tissue types
%       makec2f.m       = makes colormap for fluence rate
%
%   This example sets myname = 'skinvesseloriginal'.
%
% 7/feb/2017, add boundaryflag (see A(10)).
% 1/june/2017 , no major changes, just clean up display outputs.
% Steven L Jacques
home; clear
format compact
commandwindow

SAVEPICSON = 1;
if SAVEPICSON
    sz = 10; fz = 7; fz2 = 5; % to use savepic.m
else
    sz = 12; fz = 9; fz2 = 7; % for screen display
end

%%%% USER CHOICES <---------- you must specify -----
myname = 'skinvesseloriginal'; nm = 532;
%%%%


disp(sprintf('------ mcxyz %s -------',myname))

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
Nt = A(22);
j = 22;
for i=1:Nt
    j=j+1;
    muav(i,1) = A(j);
    j=j+1;
    musv(i,1) = A(j);
    j=j+1;
    gv(i,1) = A(j);
end

reportHmci(myname)

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

%%
x = ([1:Nx]-Nx/2-1/2)*dx;
y = ([1:Ny]-Ny/2-1/2)*dx;
z = ([1:Nz]-1/2)*dz;
ux = [2:Nx-1];
uy = [2:Ny-1];
uz = [2:Nz-1];
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;
xmin = min(x);
xmax = max(x);
xdiff = xmax-xmin;

%% Look at structure, Tzx
Tzx = reshape(T(Ny/2,:,:),Nx,Nz)';
tissue = make_TissueList_original(nm);
Nt = length(tissue);

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
N = 20; % # of beam rays drawn
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
    name = sprintf('%s_tissue',myname);
    savefig(name)
end


%% Look at Fluence Fzx @ launch point
Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source

figure(2);clf
imagesc(x,z,log10(Fzx),[.5 2.8])
hold on
text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
colormap(makec2f)
axis equal image
%axis([min(x) max(x) min(z) max(z)])
text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',time_min),...
    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%s_Fzx',myname);
    savefig(name)
end

%% look Fzy
Fzy = reshape(F(:,Nx/2,:),Ny,Nz)';

iy = round((dy*Ny/2 + 0.15)/dy);
iz = round(zs/dz);
zzs  = zs;
%Fdet = mean(reshape(Fzy(iz+[-1:1],iy+[0 1]),6,1));

figure(3);clf
imagesc(y,z,log10(Fzy),[.5 2.8])
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

drawnow

disp('done')


