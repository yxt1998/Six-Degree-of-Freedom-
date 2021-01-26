function pp
clear;close;clc
time=0:0.05:5;
m=10;
r=1;
I=0.5*m*r^2*eye(3);
u=rand(4,1);
e=u/norm(u);
et=e';
w=[0;0;5];
E=[-et(2) -et(3) -et(4) et(1);et(1) -et(4) et(3) et(2);et(4) et(1) -et(2) et(3);-et(3) et(2) et(1) et(4)];
%check
A=E*E'

nt=[eye(3);zeros(1,3)];
%check
B=nt'*nt

ed=0.5*E*nt*w

%check
ww=2*nt'*E'*ed

xx0=[e(1);e(2);e(3);e(4);ed(1);ed(2);ed(3);ed(4)];
[Tout, yot]=ode45(@dil,time,xx0)
for l=1:length(Tout)
    V=[yot(l,1);yot(l,2);yot(l,3);yot(l,4)];
    vv=norm(V)
    plot(Tout,vv,'.');
end 
end
function xd=dil(time,X)
xd=zeros(8,1);
m=10;
r=1;
I=0.5*m*r^2*eye(3);

nt=[eye(3);zeros(1,3)];

z=[X(1);X(2);X(3);X(4)];

E=[-z(2) -z(3) -z(4) z(1);z(1) -z(4) z(3) z(2);z(4) z(1) -z(2) z(3);-z(3) z(2) z(1) z(4)];
L=2*nt'*E';
zd=[X(5);X(6);X(7);X(8)];
K=cross((L*zd),I*L*zd);
M=[1;2;10];
B=[I*L;z'];
C=[(M-K);(-zd'*zd)];
edd=B\C;

xd=[zd;edd];
end 