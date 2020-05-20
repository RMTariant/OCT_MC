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

%disp(sprintf('------ mcxyz %s -------',myname))

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
%zsurf = A(22);
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
    DetID = fread(fid, 'float');
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

%% Construct the signal
k=18796.99;% wavenumber for 532 nm (simulated wavelength)
j=2;%factor 2 to denote the signal forward and backward scattering
s = exp(1j.*k.*S').*sqrt(W.*L)';
s = s/length(ix);
I = (abs(0.01*s+1)).^2 - (abs(0.01*s -1) ).^2;
y = I;

%% Apply the window
M=length(y); % length of the hamming window from y
y = y.*hamming(M);
%% Zero padding and FFT
M10 = length(y).*10;
Y10 = fft(y,M10);
Fs=M10;% Sampling frequency - legth of the M10 matrix
z10 = ((0:(M10-1)).*Fs)/M10;
z10 = z10/2;
z10 = z10(1:floor(length(z10))/2);
Y10 = Y10(1:floor(length(Y10))/2);
out = db(Y10);
plot(out)
%% KE: Listing A.2 from Zhao's Code: Load the saved photons

M = length(DetL);
ix = find(DetID);

S = DetS(ix)';
W = DetW(ix)';
L = DetL(ix)';
s = exp(1j*nm*S)* sqrt(W.*L)';
s = s/length(ix);
I = (abs(0.01*s+1)).^2 - (abs(0.01*s-1)).^2;
y = I;

%% KE: Listing A.2 from Zhao's Code: Apply the window
y = y.*hamming(M);
%wvtool(hamming(800));
%hold on 

%% KE: Listing A.2 from Zhao's Code: Zero padding and FFT
M10 = length(y) * 10;
Y10 = fft(y, M10);
Fs=1000;
z10 = (0:(M10-1))*Fs/M10;
z10 = z10/2;
z10 = z10(1:floor(length(z10)/2));
Y10 = Y10(1:floor(length(Y10)/2));
out = db(Y10);
wvtool(out);
%% KE: Listing A.2 from Zhao's Code: Save the A-Scan
%OCT(:,ID) = out;
%}

