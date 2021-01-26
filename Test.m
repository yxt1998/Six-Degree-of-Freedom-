function project
clear;close;clc
time=0:0.1:5;
m=10;                           %mass of the body
r=1;                            %radius 
I=0.5*m*r^2*eye(3,3)     
u=rand(4,1);
e=u/norm(u);
w=[0;0;2];
E=[-e(2) -e(3) -e(4) e(1);e(1) -e(4) e(3) e(2);e(4) e(1) -e(2) e(3);-e(3) e(2) e(1) e(4)];
nt=[eye(3);zeros(1,3)];
ed=0.5*E*nt*w;
x0=[0 ; 2; 3; 0; 0; 0;e(1);e(2);e(3);e(4);ed(1);ed(2);ed(3);ed(4)];     %initial conidtions 
[Tout, yot]=ode45(@dill ,time,x0)
%constraints of the eulers parameter
%check
A=E*E'

%check
B=nt'*nt

ed=0.5*E*nt*w;

%check
ww=2*nt'*E'*ed;
figure(1);
hold on
for i=1:length(Tout)
  pause(0.1);
    plot3(yot(1:i,1), yot(1:i,2), yot(1:i,3),'*');  % Plotting trajectory
    title('The trajectory of the body with Cross Wind');
    view(50,50)
    grid on
    hold on;
end %plotting the trajectory
figure(2)
for q=1:length(Tout)
    ve=[yot(q,4);yot(q,5);yot(q,6)];
    VE=dot(ve,ve');
    E=0.5*m*VE+m*9.8*yot(q,3)+0.5*dot((I*w),w);
    plot(Tout,E,'.');
    title('The energy of the System');
    xlabel('Time');
    ylabel('Energy');
    fprintf('The energy of the system is:%.4f\n',E);
end %Energy of the system 
figure(3)
for l=1:length(Tout)
    Q=[yot(l,7);yot(l,8);yot(l,9);yot(l,10)];
    hold on
    plot(Tout,norm(Q),'*');
    title('The norm of the quatronions @ any given time');
    xlabel('Time');
    ylabel('Norm')
end %Norm of the system
figure(4)
for j=1:length(Tout)
epi=[yot(j,8) yot(j,9) yot(j,10)]';
esk=[0 -yot(j,10) yot(j,9);yot(j,10) 0 -yot(j,8);-yot(j,9) yot(j,8) 0];%skew matrix
RNA=(eye(3)*(1-(2*(epi'*epi))))+(2*(epi*epi'))+2*yot(j,7)*esk;
T=[RNA [yot(j,1) yot(j,2) yot(j,3)]';zeros(1,3) 1];
pna=[(1)^2 0 0 0 0;0 0 (1)^2 0 0;0 0 0 0 (1)^2;1 1 1 1 1];
P=(T*pna);
hold on
grid on
pause(0.1);
plot3(pna(1,:), pna(2,:),pna(3,:),'-k');
plot3(P(1,:),P(2,:),P(3,:),'-r');
title('Rotation Of the frame');
view(50,50)
end 
end 
function Xd=dill(time,X)
cd=0.005;
A=0.2;
p=1;
Xd=zeros(14,1);
Xd(1)=X(4);
Xd(2)=X(5);
Xd(3)=X(6);
%flowrate=[0;-100;0];
%V=[X(4); X(5); X(6)];
%f=cd*A*p*dot(V,V)*V/norm(V);

Xd(4)=0;%f(1)
Xd(5)=0;%f(2)
Xd(6)=0;           %making gravity zero
m=10;
r=1;
I=0.5*m*r^2*eye(3,3);
nt=[eye(3);zeros(1,3)];
z=[X(7);X(8);X(9);X(10)];

E=[-z(2) -z(3) -z(4) z(1);z(1) -z(4) z(3) z(2);z(4) z(1) -z(2) z(3);-z(3) z(2) z(1) z(4)];
L=2*nt'*E';
zd=[X(11);X(12);X(13);X(14)];
K=cross((L*zd),I*(L*zd));
M=[0;0;0];
B=[I*L;z'];
C=[(M-K);(-zd'*zd)];
edd=B\C;
end 