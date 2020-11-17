% RMT
% Compared path length with oct process data

%% Select the oct data and path length

z;
OCTt = OCT;
r = S2(1,:)*1.1+S2(2,:)*1.2; %add the effct of refractive index
W;
L;

%Normalize the A-line
dx = z(2)-z(1);
OCTt = OCTt./sum(OCTt*dx*2,1);
OCT = mean(OCTt,2);

%Normalize and weigth to A-line
[counts edges] = histcounts(r,7592);
for n = 1:length(counts)
    sel = and((r > edges(n)),(r < edges(n+1)));
    counts(n) = sum(W(sel).*L(sel));
end

dx2 = edges(2)-edges(1);
counts = counts./sum(counts*dx2);

plot(edges(1:end-1)+dx2/2,counts,z*2,OCT)
legend('Aline','Photon pathlength')
xlabel('Depth (mm)')
ylabel('Intensity')

figure
semilogy(edges(1:end-1)+dx2/2,counts,z*2,OCT)
legend('Aline','Photon pathlength')
xlabel('Depth (mm)')
ylabel('Intensity (dB)')