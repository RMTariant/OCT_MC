
filepath = 'C:\Users\raphi\Documents\ubuntu_share\';
cd(filepath)
myname = 'test2';

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%Number of photon loaded at the same time
nph = 100000; %%%%CHOOSE VALUE%%%%
%Limit the total number of loaded photons if you don't want to analyse all the data
%Input 0 if you want to load the whole data
maxnph = 0;
%OCT wavelengths IN CENTIMETERS
lambda_start = 1150e-7; %%%%CHOOSE VALUE%%%%
lambda_stop = 1450e-7; %%%%CHOOSE VALUE%%%%
%Number of sample point of the OCT wavelength width
samplePoints= 2048; %%%%CHOOSE VALUE%%%%
%Choose to apply electric filter. Put a high value if none
maxDepth = 0.25; %%%%CHOOSE VALUE%%%%
%Compress image for refractive index
n_cor = 1;
%Chosse refractive of the different mediums. Index. Can be a function of the wavelength
rn = ones(samplePoints,1); %%%%CHOOSE VALUE%%%%
rn(1:samplePoints,1) = 1; %%%%CHOOSE VALUE%%%%
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

%% Load the fluence
filename = sprintf('%s_F.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)

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

%Apply low pass filter
ksampling = 2*pi/(k(1)-k(2));
rawAline = lowpass(I.*hann(length(k)),maxDepth*2*n_cor,ksampling);

%Calculate Aline
M10 = length(rawAline(:,1))*10;
OCT = abs(fft(rawAline,M10)).^2;
OCT = OCT(1:floor(length(OCT(:,1))/2)+1,:);
OCT(2:end-1) = 2*OCT(2:end-1);


%% Displaying image

z = pi/(k(1)-k(2))/(M10/length(rawAline(:,1)))*(0:(length(k)/2))/length(k)./n_cor;
x = linspace(-radius,radius,Ndetectors);

% figure
% imagesc(x,z,db(OCT))
% title('Refractive index equal to 1')
% xlabel('Position [cm]')
% ylabel('Depth [cm]')

figure
imagesc(x,z,db(OCT))
title('')
xlabel('Position [cm]')
ylabel('Depth [cm]')
caxis([-261 305])

%% Average everything!
figure
plot(z,db(mean(OCT,2)))
