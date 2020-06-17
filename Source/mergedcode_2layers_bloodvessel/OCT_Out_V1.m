%% Parameters

%File name
myname = 'skinvessel2layersb'; %%%%CHOOSE VALUE%%%%
%OCT wavelengths IN CENTIMETERS
lambda_start = 1220e-7; %%%%CHOOSE VALUE%%%%
lambda_stop = 1380e-7; %%%%CHOOSE VALUE%%%%
%Number of sample point of the OCT wavelength width
samplePoints= 512; %%%%CHOOSE VALUE%%%%

%% Load data

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);
Nx = A(2);
Ny = A(3);
Nz = A(4);
dx = A(5);
dy = A(6);
dz = A(7);
radius = A(20);

% Load path lengths of detected photons DetS
filename = sprintf('%s_DetS.bin',myname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    DetS = fread(fid, 'float');
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
W = DetW(ix)';
L = DetL(ix)';
ID = DetID(ix)'; %Raphael
Ndetector = max(ID);
%% k sampling, (needs to be corrected for proper values) Raphael

k_start = 2*pi/lambda_stop;
k_stop = 2*pi/lambda_start;
k = linspace(k_stop,k_start,samplePoints)';


%% Construct signal
%Building photon phase and amplitude
s = exp(1i.*k.*S).*sqrt(W.*L);

%Adding all the photons from a same detector together
s_sample = zeros(length(k),max(ID));
nphoton_dec = zeros(1,max(ID));
for n = 1:max(ID)
    s_sample(:,n) = sum(s(:,find(ID == n)),2);
    nphoton_dec(n) = sum(ID == n);
end

%Calculate the average amplitude of the sample arm signal to calibrate
%reference arm signal amplitude
ref_amp = max(mean(abs(s_sample),1)); %RMT: To be update. Not sure if ok.
%ref_amp = mean(abs(s_sample),'all'); alternatively

%Make the interference with the reference arm and calculate the intensity
I = (abs(s_sample + ref_amp)).^2 - (abs(s_sample - ref_amp) ).^2;

%% Processing the OCT signal
OCT = abs(fft(I.*hann(length(k)))).^2;
OCT = OCT(1:floor(length(OCT(:,1))/2)+1,:);
OCT(2:end-1) = 2*OCT(2:end-1);

%% Displaying image

z = 1/(k(1)-k(2))*(0:(length(k)/2))/length(k);
x = linspace(-radius,radius,Ndetector);

imagesc(x,z,db(OCT))
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







