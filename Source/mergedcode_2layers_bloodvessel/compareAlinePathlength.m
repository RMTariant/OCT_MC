% RMT
% Compared path length with oct process data

%% Select the oct data and path length

z;
OCTt = OCT;
r = S2*1; %add the effct of refractive index
W;
L;

%Normalize the A-line
dx = z(2)-z(1);
OCTt = OCTt./sum(OCTt*dx*2,1);
OCT = mean(OCTt,2);

%Normalize and weigth to A-line
[counts edges] = histcounts(r,1000);
for n = 1:length(counts)
    sel = and((r > edges(n)),(r < edges(n+1)));
    counts(n) = sum(W(sel).*L(sel));
end

dx = edges(2)-edges(1);
counts = counts./sum(counts*dx);

plot(z*2,OCT,edges(1:end-1)+dx/2,counts)
legend('Aline','Photon pathlength')
xlabel('Depth (mm)')
ylabel('Intensity')

figure
semilogy(z*2,OCT,'.',edges(1:end-1)+dx/2,counts,'.')
legend('Aline','Photon pathlength')
xlabel('Depth (mm)')
ylabel('Intensity (dB)')