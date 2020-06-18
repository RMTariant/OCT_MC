%% Parameters

%File name
myname = 'skinvessel2layersb'; %%%%CHOOSE VALUE%%%%
%OCT wavelengths IN CENTIMETERS
lambda_start = 1100e-7; %%%%CHOOSE VALUE%%%%
lambda_stop = 1400e-7; %%%%CHOOSE VALUE%%%%
%Number of sample point of the OCT wavelength width
samplePoints= 512; %%%%CHOOSE VALUE%%%%
%Chosse refractive of the different mediums. Index. Can be a function of the wavelength
rn = ones(samplePoints,2); %%%%CHOOSE VALUE%%%%
rn(1:samplePoints,1) = 1.3; %%%%CHOOSE VALUE%%%%
rn(1:samplePoints,2) = 1; %%%%CHOOSE VALUE%%%%

%% Load data

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);
Ndetectors = A(3);
Nx = A(4);
Ny = A(5);
Nz = A(6);
dx = A(7);
dy = A(8);
dz = A(9);
radius = A(22);
Nt = A(25);

% Load path lengths of detected photons DetS
filename = sprintf('%s_DetS.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetS = fread(fid, 'float');
    fclose(fid);
toc

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

%% Remove the outliers using .9 quantile of L 

L_threshold = quantile(DetL,0.9);
ix = find(DetL < L_threshold );

% end of removing outliers
%% Load the saved photons
S = DetS(ix)'; % row vectors
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
s = exp(1i.*k.*S).*sqrt(W.*L);
%Alternative with refractive index
kns2 = zeros(samplePoints,length(W));
for n = 1:Nt
    kns2 = kns2 + kn(:,n).*S2(n,:);
end
s2 = exp(1i.*kns2).*sqrt(W.*L);
clear kns2

%Adding all the photons from a same detector together
s_sample = zeros(samplePoints,Ndetectors);
nphoton_dec = zeros(1,Ndetectors);
for n = 1:Ndetectors
    s_sample(:,n) = sum(s(:,find(ID == n)),2);
    nphoton_dec(n) = sum(ID == n);
end
clear s
s_sample2 = zeros(samplePoints,Ndetectors);
for n = 1:Ndetectors
    s_sample2(:,n) = sum(s2(:,find(ID == n)),2);
end
clear s2

%Calculate the average amplitude of the sample arm signal to calibrate
%reference arm signal amplitude
ref_amp = max(mean(abs(s_sample),1)); %RMT: To be update. Not sure if ok.
%ref_amp = mean(abs(s_sample),'all'); alternatively

%Make the interference with the reference arm and calculate the intensity
I = (abs(s_sample + ref_amp)).^2 - (abs(s_sample - ref_amp) ).^2;
I2 = (abs(s_sample2 + ref_amp)).^2 - (abs(s_sample2 - ref_amp) ).^2;

%% Processing the OCT signal
OCT = abs(fft(I.*hann(length(k)))).^2;
OCT = OCT(1:floor(length(OCT(:,1))/2)+1,:);
OCT(2:end-1) = 2*OCT(2:end-1);

OCT2 = abs(fft(I2.*hann(length(k)))).^2;
OCT2 = OCT2(1:floor(length(OCT2(:,1))/2)+1,:);
OCT2(2:end-1) = 2*OCT2(2:end-1);

%% Displaying image

z = pi/(k(1)-k(2))*(0:(length(k)/2))/length(k);
x = linspace(-radius,radius,Ndetectors);

figure
title('Refractive index equal to 1')
imagesc(x,z,db(OCT))
xlabel('Position [cm]')
ylabel('Depth [cm]')

figure
title('Example of impact of the refractive index')
imagesc(x,z,db(OCT2))
xlabel('Position [cm]')
ylabel('Depth [cm]')

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







