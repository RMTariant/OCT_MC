%% Remove the outliers using .9 quantile of L 
L_threshold = quantile(DetL,0.9);
ix = find(DetL < L_threshold );
% end of removing outliers
%% Load the saved photons
S = DetS(ix)'; % row vectors
W = DetW(ix)';
L = DetL(ix)';
ID = DetID(ix)'; %Raphael

%% k sampling, (needs to be corrected for proper values) Raphael
lambda_start = 1220e-7;
lambda_stop = 1380e-7;
k_start = 2*pi/lambda_stop;
k_stop = 2*pi/lambda_start;
k = linspace(k_stop,k_start,1024)';

%% Construct the signal
% k=18796.99;% wavenumber for 532 nm (simulated wavelength)
% j=2;%factor 2 to denote the signal forward and backward scattering
% s = exp(1j.*k.*S').*sqrt(W.*L)';
% s = s/length(ix);
% I = (abs(0.01*s+1)).^2 - (abs(0.01*s -1) ).^2;
% y = I;

%% Construct signal
%Building photon phase and amplitude
s = exp(1i.*k.*S).*sqrt(W.*L);

s = exp(1i.*k.*S).*sqrt(W);

%Adding all the photons from a same detector together
s_sample = zeros(length(k),max(ID));
nphoton_dec = zeros(1,max(ID));
for n = 1:max(ID)
    s_sample(:,n) = sum(s(:,find(ID == n)),2);
    nphoton_dec(n) = sum(ID == n);
end

%Calculate the average amplitude of the sample arm signal to calibrate
%reference arm signal amplitude
ref_amp = max(mean(abs(s_sample),1));
%ref_amp = mean(abs(s_sample),'all'); alternatively

%Make the interference with the reference arm and calculate the intensity
I = (abs(s_sample + ref_amp)).^2 - (abs(s_sample - ref_amp) ).^2;

%% Apply the window
% M=length(y); % length of the hamming window from y
% y = y.*hamming(M);
%% Zero padding and FFT
% M10 = length(y).*10;
% Y10 = fft(y,M10);
% Fs=M10;% Sampling frequency - legth of the M10 matrix
% z10 = ((0:(M10-1)).*Fs)/M10;
% z10 = z10/2;
% z10 = z10(1:floor(length(z10))/2);
% Y10 = Y10(1:floor(length(Y10))/2);
% out = db(Y10);
% plot(out)
% % Working the code below to see the A-line

%% Processing the OCT signal
OCT = abs(fft(I.*hann(length(k)))).^2;
OCT = OCT(1:floor(length(OCT(:,1))/2)+1,:);
OCT(2:end-1) = 2*OCT(2:end-1);

imagesc(db(OCT))

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

%% Debugging
n = 144;
s_sample = sum(s(:,find(ID==n)),2);
S_sample = S(:,find(ID==n));
test = L(ID==n);
hist(test)
ref_amp = max(mean(abs(s_sample),1));

I = (abs(s_sample + ref_amp)).^2 - (abs(s_sample - ref_amp) ).^2;


%% Inspect where the photons come from.
L_threshold = quantile(DetL,0.9);
ix = find(DetL < L_threshold );
S = DetS(ix)'; % row vectors
W = DetW(ix)';
L = DetL(ix)';
ID = DetID(ix)'; %Raphael

zedges = linspace(0,0.8e-3,513);
qty = W.*L;
z = discretize(S,zedges);

I = zeros(512);
for n = 1:length(ID)
    x = ID(n);
    if qty2(n) == nan
        
    else
        I(z,x) = qty(n)+I(z,x);
    end
    n
end

I2 = zeros(512);
for nz = 1:512
    for nx = 1:512
        I2(nz,nx) = sum(qty(and(nz == z, nx == ID)));   
    end
    nz
end

%
S2 = S/2;
Qty = W.*L;
Qty = L;
Qty = W;
Qty = ones(size(L));
zedges = linspace(0,max(S2),201);
z = discretize(S2,zedges);
for n = 1:200
    test(n) = sum(Qty(z == n));
end
plot((zedges(1:end-1)+(zedges(2)-zedges(1))/2)*1e3,test)
ylabel('Intensity [a.u.]')
xlabel('Depth [mm]')

