% function to plot nodes in 3-D space
function plot_nodes(setAppropriatePositions,targetpoints,population,sensingRange,str)

%This is for ploting in the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    F9=plot3(setAppropriatePositions(:,1),setAppropriatePositions(:,2),setAppropriatePositions(:,3),'.','color','r');
    hold on
    F7=plot3(targetpoints(:,1),targetpoints(:,2),targetpoints(:,3),'*','color','b');
    hold on
     %plot the circular transmission range for population
    for ii=1:size(population,1)
        [x1,y1,z1]=sph(population(ii,1),population(ii,2),population(ii,3),sensingRange);
        surf(x1,y1,z1,'FaceAlpha',.3, 'EdgeColor', 'none');
        alpha 0.3;
        plot3(x1(z1==0),y1(z1==0),z1(z1==0),'k--','LineWidth',1)
        [x2,y2,z2]=sph2cart(linspace(0,2*pi,100),-pi/7,1);
        plot3(x2,y2,repelem(z2,100),'k--','LineWidth',1)
        hold on
    end
    axis on
    xlabel('x(m)')
    ylabel('y(m)')
    title(str)
    axis equal
end


function [x1,y1,z1]=sph(x,y,z,r)
% Set the number of points to use for plotting the sphere
n = 50;

% Generate the x, y, z coordinates for a unit sphere
[x1, y1, z1] = sphere(n);

% Scale the sphere to the desired radius
x1 = x1 * r + x ;
y1 = y1 * r + y;
z1 = z1 * r + z;

% Plot the sphere
%figure;
%surf(x1, y1, z1);
%axis equal;
end




