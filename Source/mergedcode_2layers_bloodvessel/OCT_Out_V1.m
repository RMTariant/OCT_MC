%% Remove the outliers using .9 quantile of L
automatic 
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
% Working the code below to see the A-line
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
