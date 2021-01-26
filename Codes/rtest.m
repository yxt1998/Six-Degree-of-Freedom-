clc
clear
m=4;
r=0.05;
time=0:0.1:5;
I=0.5*m*(r^2)*eye(3)
u=[1;2;3;4]
e=u/norm(u)
w=[0;0;5]
E=[-e(2) -e(3) -e(4) e(1);e(1) -e(4) e(3) e(2);e(4) e(1) -e(2) e(3);-e(3) e(2) e(1) e(4)]
nt=[eye(3);zeros(1,3)]
ed=0.5*E*nt*w
L=2*nt'*E'
K=cross((L*ed),I*(L*ed))
M=[0;0;2]
B=[I*L;e']
C=[(M-K);(-ed'*ed)]
edd=C\B