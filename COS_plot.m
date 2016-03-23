



function [ axes ] = COS_plot( wn1,wn2, COS, levels )
%COS_PLOT Summary of this function goes here
%   Detailed explanation goes here
fig = figure();
[x, y] = meshgrid(wn1,wn2);

abs_max = max(max(abs(COS)));

if length(levels) == 1 
    contour_levels = linspace(-abs_max, abs_max, 2*levels+1);
    last_neg = contour_levels(levels);
else
    contour_levels = [-levels,levels] * abs_max;
    num_levels = length(levels);
    last_neg = max(contour_levels(contour_levels<0));
end 

colormap gray;

size(y)
COS_neg = COS;

COS_neg(COS > last_neg) = 0;

hold on;

contourf(x,y,COS_neg,50,'edgecolor','none');
caxis([min(min(contour_levels)), last_neg]);



contour(x,y,COS,contour_levels(contour_levels >0), 'edgecolor','black');
contour(x,y,COS,contour_levels(contour_levels <0), 'edgecolor','black','linestyle','--');


set(gca, 'XDir', 'reverse', 'YDir','reverse');

axis image;


xlabel('\nu_1')
ylabel('\nu_2')
hold off;
end

