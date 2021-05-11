function createfigure(cdata1)
%CREATEFIGURE(cdata1)
%  CDATA1:  image cdata

%  Auto-generated by MATLAB on 09-Apr-2020 12:54:01

% Create figure
figure('PaperSize',[4 3],...
    'Colormap',[0.05 0.05 0.05;0.1 0.1 0.1;0.15 0.15 0.15;0.2 0.2 0.2;0.25 0.25 0.25;0.3 0.3 0.3;0.35 0.35 0.35;0.4 0.4 0.4;0.45 0.45 0.45;0.5 0.5 0.5;0.45 0.45 0.55;0.4 0.4 0.6;0.35 0.35 0.65;0.3 0.3 0.7;0.25 0.25 0.75;0.2 0.2 0.8;0.15 0.15 0.85;0.1 0.1 0.9;0.05 0.05 0.95;0 0 1;0 0.0833333333333333 1;0 0.166666666666667 1;0 0.25 1;0 0.333333333333333 1;0 0.416666666666667 1;0 0.5 1;0 0.583333333333333 1;0 0.666666666666667 1;0 0.75 1;0 0.833333333333333 1;0 0.916666666666667 1;0 1 1;0.0769230769230769 1 0;0.153846153846154 1 0;0.230769230769231 1 0;0.307692307692308 1 0;0.384615384615385 1 0;0.461538461538462 1 0;0.538461538461538 1 0;0.615384615384615 1 0;0.692307692307692 1 0;0.769230769230769 1 0;0.846153846153846 1 0;0.923076923076923 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 1 0;1 0.909090909090909 0;1 0.818181818181818 0;1 0.727272727272727 0;1 0.636363636363636 0;1 0.545454545454545 0;1 0.454545454545455 0;1 0.363636363636364 0;1 0.272727272727273 0;1 0.181818181818182 0;1 0.0909090909090909 0;1 0 0],...
    'OuterPosition',[-6.33333333333333 33.6666666666667 1294.66666666667 694.666666666667]);

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create image
image([-0.04975 -0.04925 -0.04875 -0.04825 -0.04775 -0.04725 -0.04675 -0.04625 -0.04575 -0.04525 -0.04475 -0.04425 -0.04375 -0.04325 -0.04275 -0.04225 -0.04175 -0.04125 -0.04075 -0.04025 -0.03975 -0.03925 -0.03875 -0.03825 -0.03775 -0.03725 -0.03675 -0.03625 -0.03575 -0.03525 -0.03475 -0.03425 -0.03375 -0.03325 -0.03275 -0.03225 -0.03175 -0.03125 -0.03075 -0.03025 -0.02975 -0.02925 -0.02875 -0.02825 -0.02775 -0.02725 -0.02675 -0.02625 -0.02575 -0.02525 -0.02475 -0.02425 -0.02375 -0.02325 -0.02275 -0.02225 -0.02175 -0.02125 -0.02075 -0.02025 -0.01975 -0.01925 -0.01875 -0.01825 -0.01775 -0.01725 -0.01675 -0.01625 -0.01575 -0.01525 -0.01475 -0.01425 -0.01375 -0.01325 -0.01275 -0.01225 -0.01175 -0.01125 -0.01075 -0.01025 -0.00975 -0.00925 -0.00875 -0.00825 -0.00775 -0.00725 -0.00675 -0.00625 -0.00575 -0.00525 -0.00475 -0.00425 -0.00375 -0.00325 -0.00275 -0.00225 -0.00175 -0.00125 -0.00075 -0.00025 0.00025 0.00075 0.00125 0.00175 0.00225 0.00275 0.00325 0.00375 0.00425 0.00475 0.00525 0.00575 0.00625 0.00675 0.00725 0.00775 0.00825 0.00875 0.00925 0.00975 0.01025 0.01075 0.01125 0.01175 0.01225 0.01275 0.01325 0.01375 0.01425 0.01475 0.01525 0.01575 0.01625 0.01675 0.01725 0.01775 0.01825 0.01875 0.01925 0.01975 0.02025 0.02075 0.02125 0.02175 0.02225 0.02275 0.02325 0.02375 0.02425 0.02475 0.02525 0.02575 0.02625 0.02675 0.02725 0.02775 0.02825 0.02875 0.02925 0.02975 0.03025 0.03075 0.03125 0.03175 0.03225 0.03275 0.03325 0.03375 0.03425 0.03475 0.03525 0.03575 0.03625 0.03675 0.03725 0.03775 0.03825 0.03875 0.03925 0.03975 0.04025 0.04075 0.04125 0.04175 0.04225 0.04275 0.04325 0.04375 0.04425 0.04475 0.04525 0.04575 0.04625 0.04675 0.04725 0.04775 0.04825 0.04875 0.04925 0.04975],...
    [0.00025 0.00075 0.00125 0.00175 0.00225 0.00275 0.00325 0.00375 0.00425 0.00475 0.00525 0.00575 0.00625 0.00675 0.00725 0.00775 0.00825 0.00875 0.00925 0.00975 0.01025 0.01075 0.01125 0.01175 0.01225 0.01275 0.01325 0.01375 0.01425 0.01475 0.01525 0.01575 0.01625 0.01675 0.01725 0.01775 0.01825 0.01875 0.01925 0.01975 0.02025 0.02075 0.02125 0.02175 0.02225 0.02275 0.02325 0.02375 0.02425 0.02475 0.02525 0.02575 0.02625 0.02675 0.02725 0.02775 0.02825 0.02875 0.02925 0.02975 0.03025 0.03075 0.03125 0.03175 0.03225 0.03275 0.03325 0.03375 0.03425 0.03475 0.03525 0.03575 0.03625 0.03675 0.03725 0.03775 0.03825 0.03875 0.03925 0.03975 0.04025 0.04075 0.04125 0.04175 0.04225 0.04275 0.04325 0.04375 0.04425 0.04475 0.04525 0.04575 0.04625 0.04675 0.04725 0.04775 0.04825 0.04875 0.04925 0.04975 0.05025 0.05075 0.05125 0.05175 0.05225 0.05275 0.05325 0.05375 0.05425 0.05475 0.05525 0.05575 0.05625 0.05675 0.05725 0.05775 0.05825 0.05875 0.05925 0.05975 0.06025 0.06075 0.06125 0.06175 0.06225 0.06275 0.06325 0.06375 0.06425 0.06475 0.06525 0.06575 0.06625 0.06675 0.06725 0.06775 0.06825 0.06875 0.06925 0.06975 0.07025 0.07075 0.07125 0.07175 0.07225 0.07275 0.07325 0.07375 0.07425 0.07475 0.07525 0.07575 0.07625 0.07675 0.07725 0.07775 0.07825 0.07875 0.07925 0.07975 0.08025 0.08075 0.08125 0.08175 0.08225 0.08275 0.08325 0.08375 0.08425 0.08475 0.08525 0.08575 0.08625 0.08675 0.08725 0.08775 0.08825 0.08875 0.08925 0.08975 0.09025 0.09075 0.09125 0.09175 0.09225 0.09275 0.09325 0.09375 0.09425 0.09475 0.09525 0.09575 0.09625 0.09675 0.09725 0.09775 0.09825 0.09875 0.09925 0.09975],...
    cdata1,'CDataMapping','scaled');

% Create text
text('FontSize',7,'String','log_{10}( A )','Position',[0.0597 -0.00374 0]);

% Create text
text('FontSize',5,'String','runtime = 10.0 min',...
    'Position',[-0.0597 -0.00773 0]);

% Create ylabel
ylabel('z [cm]');

% Create xlabel
xlabel('x [cm]');

% Create title
title('Deposition A [W/cm^3/W.delivered] ','FontSize',7);

box(axes1,'on');
axis(axes1,'tight');
axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1],'Layer','top','PlotBoxAspectRatio',...
    [1 1 20]);
% Create colorbar
colorbar('peer',axes1);
