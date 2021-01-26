clc
clear
r = 1;
[x,y,z] = sphere(50);
x0 = 0; y0 = 0; z0 = 0;
x = x*r + x0;
y = y*r + y0;
z = z*r + z0;
figure
lightGrey =0.8* [0 0.5 0.5]; % It looks better if the lines are lighter
surface(x,y,z,'FaceColor', 'none','EdgeColor',lightGrey)
hold on
grid on