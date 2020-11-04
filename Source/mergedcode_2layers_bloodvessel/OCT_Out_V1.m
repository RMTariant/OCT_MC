%% Parameters

%File name
myname = 'attenuationtest'; %%%%CHOOSE VALUE%%%%
%Number of photon loaded at the same time
nph = 40000; %%%%CHOOSE VALUE%%%%
%Limit the total number of loaded photons
maxnph = 1000000;
%OCT wavelengths IN CENTIMETERS
lambda_start = 900e-7; %%%%CHOOSE VALUE%%%%
lambda_stop = 1600e-7; %%%%CHOOSE VALUE%%%%
%Number of sample point of the OCT wavelength width
samplePoints= 2048; %%%%CHOOSE VALUE%%%%
%Choose to apply electric filter. Put a high value if none
maxDepth = 1.12; %%%%CHOOSE VALUE%%%%
%Remove the photon with very high likelihood
L_filter = 1; %%%%CHOOSE VALUE%%%%
%Compress image for refractive index
n_cor = 1;
%Chosse refractive of the different mediums. Index. Can be a function of the wavelength
rn = ones(samplePoints,5); %%%%CHOOSE VALUE%%%%
rn(1:samplePoints,1) = 1.0; %%%%CHOOSE VALUE%%%%
rn(1:samplePoints,2) = 1.0; %%%%CHOOSE VALUE%%%%
rn(1:samplePoints,3) = 1.0; %%%%CHOOSE VALUE%%%%
rn(1:samplePoints,4) = 1.0; %%%%CHOOSE VALUE%%%%
noise_amp = 0; %%%Choose amplitude of noise%%%



%% Load data

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);
Ndetectors = A(4);
Nx = A(5);
Ny = A(6);
Nz = A(7);
dx = A(8);
dy = A(9);
dz = A(10);
radius = A(23);
Nt = A(26);

% Load path lengths of detected photons DetS
% filename = sprintf('%s_DetS.bin',myname);
% disp(['loading ' filename])
% tic
%     fid = fopen(filename, 'rb');
%     DetS = fread(fid, 'float');
%     fclose(fid);
% toc

% Load path lengths of detected photons DetS
filename = sprintf('%s_DetS2.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetS2 = fread(fid, 'float');
    fclose(fid);
toc

% Load weights of detected photons DetW
filename = sprintf('%s_DetW.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetW = fread(fid, 'float');
    fclose(fid);
toc

% Load likelihood ratios of detected photons DetL
filename = sprintf('%s_DetL.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetL = fread(fid, 'float');
    fclose(fid);
toc

% Load IDs of detected photons DetID
filename = sprintf('%s_DetID.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetID = fread(fid, 'int');
    fclose(fid);
toc

% Load E probing varaible
filename = sprintf('%s_DetE.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetE = fread(fid, 'float');
    fclose(fid);
toc

%Limit the number of photons
maxnph = min([length(DetL) maxnph]);
DetS2 = DetS2(1:maxnph);
DetW = DetW(1:maxnph);
DetL = DetL(1:maxnph);
DetID = DetID(1:maxnph);
DetE = DetE(1:maxnph);

%% Remove the outliers using .9 quantile of L 

L_threshold = quantile(DetL,L_filter);
ix = find(DetL < L_threshold );
% end of removing outliers
%% Load the saved photons
% S = DetS(ix)'; % row vectors
S2 = reshape(DetS2,[Nt,length(DetL)]);
S2 = S2(:,ix);
W = DetW(ix)';
L = DetL(ix)';
ID = DetID(ix)'; %Raphael
%Ndetectors = max(ID);
%% k sampling, (needs to be corrected for proper values) Raphael

k_start = 2*pi/lambda_stop;
k_stop = 2*pi/lambda_start;
k = linspace(k_stop,k_start,samplePoints)';
kn = ones(samplePoints,Nt);
for n = 1:Nt
    kn(:,n) = k.*rn(:,n);
end

%% Construct signal
%Building photon phase and amplitude
% s = exp(1i.*k.*S).*sqrt(W.*L);
%Alternative with refractive index
%Adding all the photons from a same detector together
% s_sample = zeros(samplePoints,Ndetectors);
% nphoton_dec = zeros(1,Ndetectors);
% for n = 1:Ndetectors
%     s_sample(:,n) = sum(s(:,find(ID == n)),2);
%     nphoton_dec(n) = sum(ID == n);
% end
% clear s

sig = zeros(samplePoints,Ndetectors);
if length(W)> nph
    p = floor(length(W)/nph);
else
    p = 1;
    nph = length(W);
end
qn = length(W) - nph*p;
if qn >= 0
    q = 1;
end

for m = 1:p
    sel = (m-1)*nph+1:m*nph;
    kns2 = zeros(samplePoints,nph);
    for n = 1:Nt
        kns2 = kns2 + kn(:,n).*S2(n,sel);
    end
    
    s2 = exp(1i.*kns2).*sqrt(W(sel).*L(sel));
    clear kns2
    
    s_sample = zeros(samplePoints,Ndetectors);
    for n = 1:Ndetectors
        s_sample(:,n) = sum(s2(:,find(ID(sel) == n)),2);
    end
    clear s2
    sig = sig + s_sample;
    clear s_sample
end
if q == 1
    sel = m*nph+1:length(W);
    kns2 = zeros(samplePoints,qn);
    for n = 1:Nt
        kns2 = kns2 + kn(:,n).*S2(n,sel);
    end
    
    s2 = exp(1i.*kns2).*sqrt(W(sel).*L(sel));
    clear kns2
    
    s_sample = zeros(samplePoints,Ndetectors);
    for n = 1:Ndetectors
        s_sample(:,n) = sum(s2(:,find(ID(sel) == n)),2);
    end
    clear s2
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
% OCT = abs(fft(I.*hann(length(k)))).^2;
% OCT = OCT(1:floor(length(OCT(:,1))/2)+1,:);
% OCT(2:end-1) = 2*OCT(2:end-1);

ksampling = 2*pi/(k(1)-k(2));
rawAline = lowpass(I.*hann(length(k)),maxDepth*2*n_cor,ksampling);
OCT = abs(fft(rawAline)).^2;
%OCT = abs(fft(I.*hann(length(k)))).^2;

OCT = OCT(1:floor(length(OCT(:,1))/2)+1,:);
OCT(2:end-1) = 2*OCT(2:end-1);

%% Displaying image

z = pi/(k(1)-k(2))*(0:(length(k)/2))/length(k)./n_cor;
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

%% Save the A- scan
% ID=DetID;
% OCT(:,ID) = out;
% %% Transformation of the output - Page no:73; section:6.3
% Y_Len = zeros(size(out,1),N);
% Y_Type = zeros(size(out,1),N);
%         for ID = 1:size(out,1)
%             y = OCTY_all(out,:);
%             b = diff(y);
%             c = find(b);
%             d = diff([0,c,length(y]);
%             Y_Len(ID,1:length(d))=d;
%             Y_Type(ID,1:length(d))=y( cumsum (d));
%       end


%% Debug







