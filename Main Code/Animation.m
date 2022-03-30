%function [outputArg1,outputArg2] = Animation(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t=130;
N = 350;
x = linspace(1,350,350);
y = zeros(1,350);
x_Bind = linspace(10,350,35);
Y_Bind = 0.0004*ones(1,35);
plot(x,y,'blue','LineWidth',20)
hold on
%plot(t, 0.0010,'.b', 'MarkerSize',140)

rectangle('Position',[t-6.5 0.0003 13 0.001],'Curvature',[1,1],'FaceColor',[0 1 0]*abs(1-mod(t,(10))/5) + [1 0 0]*abs(1-mod(t+5,(10))/5))
scatter(x_Bind,Y_Bind,95,'^','black','filled')

set(gca,'ytick',[])
ylim ([-0.004 0.01])
xlim([-0.1*N 1.1*N])

r=1; %r=0 left sided, r=1 right sided
rectangle('Position',[t-6.5+r*9 0.0005 4 0.0008],'Curvature',[1,1],'FaceColor', 'y')
rectangle('Position',[t-6.5 0.0012 14 0.0002],'Curvature',[1,1],'FaceColor', 'y')






%end

