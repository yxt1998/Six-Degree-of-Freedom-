function project
clear;close;clc
time=0:0.1:5;
m=10;                           %mass of the body
r=1;                            %radius 
I=0.5*m*r^2*eye(3,3);           %inertia of the body 
u=[1;1;1;1];                    %initial position of the axis;aligns with the reference frame
%if the body is in random postion
%u=rand(4,1);

e=u/norm(u);                    %Normalizing the Euler's parameter 
w=[1;1;1];                      %Angular Velocity

E=[-e(2) -e(3) -e(4) e(1);e(1) -e(4) e(3) e(2);e(4) e(1) -e(2) e(3);-e(3) e(2) e(1) e(4)];

nt=[eye(3);zeros(1,3)];

ed=0.5*E*nt*w;                                                          %Initial derivative of the quatronion 
x0=[0 ; 0; 0; 3; 5; 50;e(1);e(2);e(3);e(4);ed(1);ed(2);ed(3);ed(4)];     %initial conidtions 
[Tout, yot]=ode45(@dill ,time,x0);                                      %Solving using ODE$%
%constraints of the eulers parameter
%check
A=E*E';          %must be I

%check
B=nt'*nt ;       %must be I


%check
ww=2*nt'*E'*ed;     %must be input W
 figure(1)
for j=1:length(Tout)
   
pause(0.05);
epi=[yot(j,8) yot(j,9) yot(j,10)]';
esk=[0 -yot(j,10) yot(j,9);yot(j,10) 0 -yot(j,8);-yot(j,9) yot(j,8) 0];%skew matrix
RNA=(eye(3)*(1-(2*(epi'*epi))))+(2*(epi*epi'))+2*yot(j,7)*esk;
T=[RNA [yot(j,1) yot(j,2) yot(j,3)]';zeros(1,3) 1];
pna=[(1)^2 0 0 0 0;0 0 (1)^2 0 0;0 0 0 0 (1)^2;1 1 1 1 1];
P=(T*pna);
hold on
grid on
plot3(pna(1,:), pna(2,:),pna(3,:),'-k');
plot3(P(1,:),P(2,:),P(3,:),'-r');
title('Rotation Of the frame');
view(43,24)

%{
F(j) = getframe(gcf) ;
      drawnow
      % create the video writer with 1 fps
  writerObj = VideoWriter('myVideo.avi');
  writerObj.FrameRate = 03;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
%}
end %plotting rotation
figure(2)
for l=1:length(Tout)
    Q=[yot(l,7);yot(l,8);yot(l,9);yot(l,10)];
    plot(Tout,norm(Q),'*');
    title('The norm of the quatronions @ any given time');
    xlabel('Time');
    ylabel('Norm')
end %Norm of the systemw
    figure(3)
for i=1:length(Tout)

  pause(0.01)
    plot3(yot(1:i,1), yot(1:i,2), yot(1:i,3),'*');  
    view(43,24)
    grid on
    hold on;
%{
    
    F(i) = getframe(gcf) ;
      drawnow
      % create the video writer with 1 fps
  writerObj = VideoWriter('myVideo2.avi');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
   %}
end% Plotting trajectory
figure(4)
for q=1:length(Tout)
    pause(0.1)
    ve=[yot(q,4);yot(q,5);yot(q,6)];
    VE=dot(ve,ve');
    KE(q)=0.5*m*VE
    PE(q)=m*9.8*yot(q,3);
    hold on 
    grid on
%     E=0.5*m*VE+m*9.8*yot(q,3)+0.5*dot((I*w),w);
    plot(q,KE(q),'*','markersize',7);
    hold on
    plot(q,PE(q),'.','markersize',7);
    title('The  kinetic energy and potential energy of the System');
    legend('KE','PE')
    xlabel('Time');
    ylabel(' Kinetic Energy');
%     fprintf('The energy of the system is:%.4f\n',E);
end %Energy of the system 
end

%The Function just calcuates rotation
%{
function Xd=dill(time,X)

Xd=zeros(14,1);
       %making gravity zero
m=10;
r=1;
I=0.5*m*r^2*eye(3,3);
nt=[eye(3);zeros(1,3)];
z=[X(7);X(8);X(9);X(10)];
E=[-z(2) -z(3) -z(4) z(1);z(1) -z(4) z(3) z(2);z(4) z(1) -z(2) z(3);-z(3) z(2) z(1) z(4)];
L=2*nt'*E';
zd=[X(11);X(12);X(13);X(14)];
www=2*nt'*E'*zd;       %must be same as initial w if M=0
K=cross((L*zd),I*(L*zd));
M=[0;0;0];
B=[I*L;z'];
C=[(M-K);(-zd'*zd)];
edd=B\C;
Xd(1)=X(4);
Xd(2)=X(5);
Xd(3)=X(6);


Xd(4)=0;
Xd(5)=0;
Xd(6)=0;    
Xd(7:14)=[zd edd];
end
%}  
%This Function Calculates rotation @Crosswind
%{
function Xd=dill(time,X)
cd=0.005;                   %coeffecient of Drag
A=0.2;                      %Area
p=1;                        %Pressure
Xd=zeros(14,1);
Xd(1)=X(4);
Xd(2)=X(5);
Xd(3)=X(6);
flowrate=[0;-100;0];        %crosswind
V=flowrate-[X(4); X(5); X(6)];  %velocity of the body
f=cd*A*p*dot(V,V)*V/norm(V);    %final Velocity

Xd(4)=0+f(1);
Xd(5)=0+f(2);
Xd(6)=-9.81+f(3);           %making gravity zero
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
Xd(7:14)=[zd edd];
end
%}
%This function calculates the energy of the system
%
function Xd=dill(time,X)
cd=0.005;                   %coeffecient of Drag
A=0.2;                      %Area
p=1;                        %Pressure
Xd=zeros(14,1);
Xd(1)=X(4);
Xd(2)=X(5);
Xd(3)=X(6);

Xd(4)=0;
Xd(5)=0;
Xd(6)=-9.81;          
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
Xd(7:14)=[zd edd];
end
%}