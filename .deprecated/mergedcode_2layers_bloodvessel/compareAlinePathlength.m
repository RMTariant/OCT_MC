% RMT
% Compared path length with oct process data

%% Select the oct data and path length

z;
OCTt = OCT;
r = sum(S2*1,1); %add the effct of refractive index
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

%%
pL = S2(1,:)+S2(2,:)*1.4;
pL = pL/2;

dz = z(2) - z(1);
nbins = length(z);

counts = zeros(nbins,1);
for n = 1:nbins
    sel = and((pL > (dz*(n-1))),(pL < (dz*n)));
    counts(n) = sum(W(sel).*L(sel));
end

figure
plot(db(counts))

clear testW testL testQ testC
for n = 129:137
    sel = and((pL > (dz*(n-1))),(pL < (dz*n)));
    testW(n-128) = sum(W(sel));
    testL(n-128) = sum(L(sel));
    testQ(n-128) = sum(sel);
    testC(n-128) = sum(W(sel).*L(sel));
end

% plot(z*2,db(mean(OCT,2)))
% 
% %%Debug
% A530 = sum(DetID == 530);
% A531 = sum(DetID == 531);
% A532 = sum(DetID == 532);