function project
clear;close;clc
time=0:0.1:5;

x0=[0 ; 0; 0; 3; 5; 50];
[Tout, yot]=ode45(@dill ,time,x0)
figure();
hold on 
for i=1:length(Tout)
  pause(0.1);
    plot3(yot(1:i,1), yot(1:i,2), yot(1:i,3),'*');  % Plotting trajectory
    view(50,50)
    grid on
    hold on;
   
end
end

function Xd=dill(time,X)
cd=0.005;
A=0.2;
p=1;
Xd=zeros(6,1);
Xd(1)=X(4);
Xd(2)=X(5);
Xd(3)=X(6);
flowrate=[0;-100;0];
V=flowrate-[X(4); X(5); X(6)];
f=cd*A*p*dot(V,V)*V/norm(V);

Xd(4)=0+f(1);
Xd(5)=0+f(2);
Xd(6)=-9.81+f(3);
end 