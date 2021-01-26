function rotation
m=4;
r=0.05;
time=0:0.1:5;
I=0.5*m*r^2*eye(3);
u=rand(4,1);
e=u/norm(u);
w=[0;0;5];
E=[-e(2) -e(3) -e(4) e(1);e(1) -e(4) e(3) e(2);e(4) e(1) -e(2) e(3);-e(3) e(2) e(1) e(4)];
nt=[eye(3);zeros(1,3)];
ed=0.5*E*nt*w;
CC=ed'*e;
Y0=[e(1);e(2);e(3);e(4);ed(1);ed(2);ed(3);ed(4)];
[t,o]=ode45(@qa,time,Y0);
end 
function yd=qa(time,Y)
yd=zeros(8,1);
m=4;
r=0.05;
I=0.5*m*r^2*eye(3);
e=[Y(1);Y(2);Y(3);Y(4)];
E=[-e(2) -e(3) -e(4) e(1);e(1) -e(4) e(3) e(2);e(4) e(1) -e(2) e(3);-e(3) e(2) e(1) e(4)];
L=2*nt'*E';
A=[Y(5);Y(6);Y(7);yd(8)];
K=cross((L.A),I*(L.A));
M=[0;0;2];
M=(I*L)/e';
N=(M-K)/(-A'*A);
edd=inv(M)*N;
yd(1)=Y(5);
yd(2)=Y(6);
yd(3)=Y(7);
yd(4)=Y(8);
yd(5)=edd(1);
yd(6)=edd(2);
yd(7)=edd(3);
yd(8)=edd(4);
ed=[Y(5);Y(6);Y(7);Y(8)];
yd=[ed;edd];

end
