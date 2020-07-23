function tissue = make_TissueList_mergedcode_2layers_bloodvessel(nm)
%   Returns the tissue optical properties at the wavelength nm:
%       tissue = [mua; mus; g]';
%       global tissuenames(i).s
%   Uses 
%       SpectralLIB.mat

%% Load spectral library
load spectralLIB.mat
%   muadeoxy      701x1              5608  double              
%   muamel        701x1              5608  double              
%   muaoxy        701x1              5608  double              
%   muawater      701x1              5608  double              
%   musp          701x1              5608  double              
%   nmLIB         701x1              5608  double   
nm=532;
MU(:,1) = interp1(nmLIB,muaoxy,nm);
MU(:,2) = interp1(nmLIB,muadeoxy,nm);
MU(:,3) = interp1(nmLIB,muawater,nm);
MU(:,4) = interp1(nmLIB,muamel,nm);
LOADED = 1;

%% Create tissueList

j=1;
tissue(j).name  = 'blood';
tissue(j).mua = 5; %[cm^-1]
tissue(j).mus = 650; %[cm^-1]
tissue(j).g   = 0.9888;
tissue(j).n   = 1; %Refractive index
tissue(j).mirror   = 0; %Mirror reflection coeficient

j=2;
tissue(j).name  = 'standard tissue';
tissue(j).mua   = 1; %[cm^-1]
tissue(j).mus   = 10; %[cm^-1]
tissue(j).g     = 0.7;
tissue(j).n   = 1; %Refractive index
tissue(j).mirror   = 0; %Mirror reflection coeficient

disp(sprintf('---- tissueList ------ \tmua   \tmus  \tg  \tmusp'))
for i=1:length(tissue)
    disp(sprintf('%d\t%15s\t%0.4f\t%0.1f\t%0.3f\t%0.1f',...
        i,tissue(i).name, tissue(i).mua,tissue(i).mus,tissue(i).g,...
        tissue(i).mus*(1-tissue(i).g)))
end
disp(' ')

