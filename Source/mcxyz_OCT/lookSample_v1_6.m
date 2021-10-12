myname = 'Gauss_Short_5A';

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%Number of photon loaded at the same time
nph = 200000; %%%%CHOOSE VALUE%%%%
%Limit the total number of loaded photons if you don't want to analyse all the data
%Input 0 if you want to load the whole data
maxnph = 0;
%OCT wavelengths IN CENTIMETERS
lambda_start = 1150e-7; %%%%CHOOSE VALUE%%%%
lambda_stop = 1450e-7; %%%%CHOOSE VALUE%%%%
%Number of sample point of the OCT wavelength width
samplePoints= 2048; %%%%CHOOSE VALUE%%%%
%Choose to apply electric filter. Put a high value if none
maxDepth = 0; %%%%CHOOSE VALUE%%%% 0 if none
%Compress image for refractive index
n_cor = 1;
%Chosse refractive of the different mediums. Index. Can be a function of the wavelength
rn = ones(samplePoints,1); %%%%CHOOSE VALUE%%%%
rn(1:samplePoints,1) = 1.3844; %%%%CHOOSE VALUE%%%%
%rn(1:samplePoints,2) = 1.4; %%%%CHOOSE VALUE%%%%
%rn(1:samplePoints,3) = 1; %%%%CHOOSE VALUE%%%%
%rn(1:samplePoints,4) = 1; %%%%CHOOSE VALUE%%%%

%Choose amplitude of noise
noise_amp = 0; %%%%CHOOSE VALUE%%%%

%% parameters
n = 1;
time_min = A(n); n = n + 1;
a_coef = A(n); n = n + 1;
p = A(n); n = n + 1;
Ndetectors = A(n); n = n + 1;
det_radius = A(n); n = n + 1;
cos_accept = A(n); n = n + 1;
lambda = A(n); n = n + 1;
f = A(n); n = n + 1;
D = A(n); n = n + 1;
z_f_img = A(n); n = n + 1;
h_step = A(n); n = n + 1;
Nx = A(n); n = n + 1;
Ny = A(n); n = n + 1;
Nz = A(n); n = n + 1;
dx = A(n); n = n + 1;
dy = A(n); n = n + 1;
dz = A(n); n = n + 1;
mcflag = A(n); n = n + 1;
launchflag = A(n); n = n + 1;
boundaryflag = A(n); n = n + 1;
xs = A(n); n = n + 1;
ys = A(n); n = n + 1;
zs = A(n); n = n + 1;
xfocus = A(n); n = n + 1;
yfocus = A(n); n = n + 1;
zfocus = A(n); n = n + 1;
ux0 = A(n); n = n + 1;
uy0 = A(n); n = n + 1;
uz0 = A(n); n = n + 1;
radius = A(n); n = n + 1;
waist = A(n); n = n + 1;
zsurf = A(n); n = n + 1;
Nt = A(n);  n = n + 1;
j = n;
for i=1:Nt
    muav(i,1) = A(j);
    j=j+1;
    musv(i,1) = A(j);
    j=j+1;
    gv(i,1) = A(j);
    j=j+1;
    nrv(i,1) = A(j);
    j=j+1;
end

rn(1:samplePoints,1) = nrv; %%%%CHOOSE VALUE%%%%

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
%% Load ID of detected photons DetID
filename = sprintf('%s_DetID.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetID = fread(fid, 'int');
    fclose(fid);
toc
%{
%% Load the fluence
filename = sprintf('%s_F.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)
%}

%% Load the saved photons
S = DetS(:)';
W = DetW(:)';
L = DetL(:)';
ID = DetID(:)';
if maxnph ~= 0
    maxnph = min([maxnph length(W(:))]);
    S = S(1:maxnph);
    W = W(1:maxnph);
    L = L(1:maxnph);
    ID = ID(1:maxnph);
end
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

%{
%% Look at Fluence Fzx @ launch point

Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source
figure(1);
imagesc(x,z,log10(Fzx))%,[.5 2.8])
hold on
%text(max(x)*1.2,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',9)
colorbar
%set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('Fluence \phi [W/cm^2/W.delivered]')% ','fontweight','normal','fontsize',9)
txt = ['Focal length = ' num2str(f) 'cm, Beam diameter = ' num2str(D) 'cm, imaging lens position = ' num2str(z_f_img) 'cm, n = ' num2str(nrv)];
subtitle(txt)
%colormap(makec2f)
axis equal image
%axis([xmin xmax zmin zmax])
%axis([-1 1 0 5])


%% look Fzy
%Fzy = reshape(F(:,Nx/2,:),Ny,Nz)';
%% Look at Fluence Fxy @ z= 0.02
Fxy = reshape(F(:,:,Nz/2),Ny,Nx)'; % in z,x plane through source
figure(2)
imagesc(y,x,log10(Fxy))%,[.5 2.8])
hold on
colorbar
xlabel('x [cm]')
ylabel('y [cm]')
title('Fluence \phi [W/cm^2/W.delivered]')%
txt = ['Focal length = ' num2str(f) 'cm, Beam diameter = ' num2str(D) 'cm, imaging lens position = ' num2str(z_f_img) 'cm, n = ' num2str(nrv)];
subtitle(txt)
axis equal image
%axis([-0.2 0.2 -0.2 0.2])
%text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',time_min),...
    %'fontsize',fz2)
%{
if SAVEPICSON
    name = sprintf('%s_Fzx',myname);
    savefig(name)
end
    %}
%}
%% Inspect data
% int = W.*L;
% binlimits = linspace(0,dz*Nz,5000)';
% binqty = zeros(length(binlimits)-1,1);
% for n = 1:(length(binlimits)-1)
%     binqty(n) = sum((binlimits(n)<S & binlimits(n+1)>S).*int);
% end
% figure
% plot((binlimits(1:end-1)+binlimits(2:end))/2,binqty)

%% k sampling, (needs to be corrected for proper values) Raphael

k_start = 2*pi/lambda_stop;
k_stop = 2*pi/lambda_start;
k = linspace(k_stop,k_start,samplePoints)';
kn = k.*rn;
%kn = ones(samplePoints,Nt);
%for n = 1:Nt
%    kn(:,n) = k.*rn(:,n);
%end

%% Construct signal

% Load the signal in p number of groups
sig = zeros(samplePoints,Ndetectors);
if length(W)> nph
    p = ceil(length(W)/nph);
else
    p = 1;
    nph = length(W);
end

% Sum the signal one group at the time
for m = 1:p
    if m == p
        sel = (m-1)*nph+1:length(W(((m-1)*nph+1):end));
    else
        sel = (m-1)*nph+1:m*nph;
    end
    
    signal = exp(1i.*kn.*S(sel)).*sqrt(W(sel).*L(sel));
    s_sample = zeros(samplePoints,Ndetectors);
    for n = 1:Ndetectors
        s_sample(:,n) = sum(signal(:,find(ID(sel) == n)),2);
    end
    clear signal
    sig = sig + s_sample;
    clear s_sample
    m/p;
end

%Adding noise
noise1 = raylrnd(noise_amp/1.253*ones(size(sig))) .* exp(1i.*2.*pi.*rand(size(sig)));
noise2 = raylrnd(noise_amp/1.253*ones(size(sig))) .* exp(1i.*2.*pi.*rand(size(sig)));

%Calculate the average amplitude of the sample arm signal to calibrate
%reference arm signal amplitude
ref_amp = max(mean(abs(sig),1)); %RMT: To be update. Not sure if ok.
%ref_amp = mean(abs(s_sample),'all'); alternatively

%Make the interference with the reference arm and calculate the intensity
% I = (abs(s_sample + ref_amp)).^2 - (abs(s_sample - ref_amp) ).^2;
I = (abs(sig + ref_amp + noise1)).^2 - (abs(sig - ref_amp + noise2) ).^2;

%% Processing the OCT signal

%Apply low pass filter and hanning window
ksampling = 2*pi*n_cor/(k(1)-k(2));
if maxDepth == 0
    rawAline = I.*hann(length(k));
else
    rawAline = lowpass(I.*hann(length(k)),maxDepth,ksampling);
end

%Calculate Aline
M10 = length(rawAline(:,1))*10; %Zero padding
OCT = abs(fft(rawAline,M10)).^2;
OCT = OCT(1:floor(length(OCT(:,1))/2)+1,:);
OCT(2:end-1) = 2*OCT(2:end-1);


%% Displaying image
z = (0:M10/2)/M10.*ksampling;
x = linspace(-radius,radius,Ndetectors);
figure();
imagesc(x,z,db(OCT))
title('')
xlabel('Position [cm]')
ylabel('Depth [cm]')
% caxis([-261 305])

%% Average all Alines in one
z = z/100;
figure();
plot(z,db(mean(OCT,2)))
toc
D5g05Long = [z', mean(OCT,2)];