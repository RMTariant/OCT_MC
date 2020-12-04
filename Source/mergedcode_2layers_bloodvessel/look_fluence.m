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

myname = 'fluencea001p02'; nm = 532;
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


% Nphoton = A(1);
% Nx = A(2);
% Ny = A(3);
% Nz = A(4);
% dx = A(5);
% dy = A(6);
% dz = A(7);
% 
% Nt = A(22);
% j = Nt;
% for i=1:Nt
%     j=j+1;
%     muav(i,1) = A(j);
%     j=j+1;
%     musv(i,1) = A(j);
%     j=j+1;
%     gv(i,1) = A(j);
% end



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
F = F/sum(F(~isnan(F(:))));
% figure
x = (0:Nx-1)*dx;
z = (0:Nz-1)*dz;
% Fzx = reshape(F(round(Ny/2),:,:),Nx,Nz)'; % in z,x plane through source
% imagesc(x,z,log10(Fzx));
% c = colorbar;
% c.Label.String = 'decibel';

%% Look at Fluence Fzx @ launch point
Fzx = reshape(F(round(Ny/2),:,:),Nx,Nz)'; % in z,x plane through source
%Fzx_norm = Fzx - min(Fzx(:));
%Fzx_norm = Fzx ./ max(Fzx(:));
figure;clf %deletes from the current figure all graphic objects
imagesc(x,z,log10(Fzx)) % specifies image location
F3 = Fzx;
% x and z specify the locations of the corner corresponding to 
% log10(Fzx)(1,1) and log10(Fzx)(m,n)
%imagesc(x,z,log10(Fzx_norm))
% displays the data in array log10(Fzx)as an image that uses the full range
% of colors in the colormap. Each element of log10(Fzx) specifies the color of one pixel of the
% image. The resulting image is an m-by-n grid of pixels with m rows and n
% columns in log10(Fzx). The row and column indices of the elements
% determine the centers of the corresponding pixels.
% hold on
% % text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',fz)
% colorbar % displays a vertical colorbar to thhe right of the current axes or chart
% % displays the current colormap and indicate the mapping of data values
% % into the colormap
% set(gca,'fontsize',sz)
% xlabel('x [cm]')
% ylabel('z [cm]')
% title('Fluence \phi [W/cm^2/W.delivered] ','fontweight','normal','fontsize',fz)
% %colormap(makec2f)
% %colormap(map) %sets the colormap for the current figure to the colormap
% %specified by map
% axis equal image
%axis([min(x) max(x) min(z) max(z)])
%text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',time_min),...
%    'fontsize',fz2)

% if SAVEPICSON
%     name = sprintf('%s_Fzx.fig',myname);
%     savefig(name)
% end

%%
%%%%%%%%%%
myname = 'fluencea001p02'; nm = 532;
cd('C:\Users\raphi\Documents\Doctorat\Uday simulation\data')

% Load header file
filename = sprintf('%s_F.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)
F1 = reshape(F(round(Ny/2),:,:),Nx,Nz)'; % in z,x plane through source

%%%%%%%%
myname = 'fluencea09p02'; nm = 532;
cd('C:\Users\raphi\Documents\Doctorat\Uday simulation\data')

% Load header file
filename = sprintf('%s_F.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)
F2 = reshape(F(round(Ny/2),:,:),Nx,Nz)'; % in z,x plane through source

%%%%%%%%
myname = 'fluencea09p05'; nm = 532;
cd('C:\Users\raphi\Documents\Doctorat\Uday simulation\data')

% Load header file
filename = sprintf('%s_F.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)
F3 = reshape(F(round(Ny/2),:,:),Nx,Nz)'; % in z,x plane through source

%%%%%%%%
myname = 'fluencea09p02_o'; nm = 532;
cd('C:\Users\raphi\Documents\Doctorat\Uday simulation\data')

% Load header file
filename = sprintf('%s_F.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)
F4 = reshape(F(round(Ny/2),:,:),Nx,Nz)'; % in z,x plane through source

%%%%%%%%
F5 = fluence.data;
F5 = reshape(F5(:,:,51),101,101);


%% Look at attenuation at center

F1 = F1./F1(10,51);
F2 = F2./F2(10,51);
F3 = F3./F3(10,51);
F4 = F4./F4(10,51);
F5 = F5./F5(10,51);

figure
imagesc(x*10,z*10,log10(F1))
xlabel('x [mm]')
ylabel('z [mm]')
title('p=0.2, a=0.001')
h = colorbar;
set(get(h,'label'),'string','Intensity (log10)');
caxis([-5.35 0.25]) 
figure
imagesc(x*10,z*10,log10(F2))
h = colorbar;
set(get(h,'label'),'string','Intensity (log10)');
xlabel('x [mm]')
ylabel('z [mm]')
title('p=0.2, a=0.9')
caxis([-5.35 0.25]) 
figure
imagesc(x*10,z*10,log10(F3))
h = colorbar;
set(get(h,'label'),'string','Intensity (log10)');
xlabel('x [mm]')
ylabel('z [mm]')
title('p=0.5, a=0.9')
caxis([-5.35 0.25]) 
figure
imagesc(x*10,z*10,log10(F4))
h = colorbar;
set(get(h,'label'),'string','Intensity (log10)');
xlabel('x [mm]')
ylabel('z [mm]')
title('original Steve-Jacques')
caxis([-5.35 0.25]) 
figure
imagesc(x*10,z*10,log10(F5))
h = colorbar;
set(get(h,'label'),'string','Intensity (log10)');
xlabel('x [mm]')
ylabel('z [mm]')
title('MTX simulator')
caxis([-5.35 0.25]) 


figure
n = 51;
plot(x*10,log10(F1(:,n)),x*10,log10(F2(:,n)),x*10,log10(F3(:,n)),x*10,log10(F4(:,n)),x*10,log10(F5(:,n)))
legend('p=0.2, a=0.001','p=0.2, a=0.9','p=0.5, a=0.9','original Steve-Jacques','MTX simulator')
ylabel('Intensity (log10)')
xlabel('z [mm]')
title('Intensity in the middle vs depth')

figure
n = 65;
plot(x*10,log10(F1(:,n)),x*10,log10(F2(:,n)),x*10,log10(F3(:,n)),x*10,log10(F4(:,n)),x*10,log10(F5(:,n)))
legend('p=0.2, a=0.001','p=0.2, a=0.9','p=0.5, a=0.9','original Steve-Jacques','MTX simulator')
ylabel('Intensity (log10)')
xlabel('z [mm]')
title('Intensity on the side vs depth')

figure
n = 20;
plot(x*10,log10(F1(n,:)),x*10,log10(F2(n,:)),x*10,log10(F3(n,:)),x*10,log10(F4(n,:)),x*10,log10(F5(n,:)))
legend('p=0.2, a=0.001','p=0.2, a=0.9','p=0.5, a=0.9','original Steve-Jacques','MTX simulator')
ylabel('Intensity (log10)')
xlabel('x [mm]')
title('Intensity at z = 0.02 mm')

figure
n = 80;
plot(x*10,log10(F1(n,:)),x*10,log10(F2(n,:)),x*10,log10(F3(n,:)),x*10,log10(F4(n,:)),x*10,log10(F5(n,:)))
legend('p=0.2, a=0.001','p=0.2, a=0.9','p=0.5, a=0.9','original Steve-Jacques','MTX simulator')
ylabel('Intensity (log10)')
xlabel('x [mm]')
title('Intensity at z = 0.08 mm')



