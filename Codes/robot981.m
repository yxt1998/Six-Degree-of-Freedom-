function robot981
SolveOrdinaryDifferentialEquations
% File  robot981.m  created by Autolev 4.0 on Wed Aug 14 00:31:43 2019


%===========================================================================
function VAR = ReadUserInput
global   LA LB LC LD MB MC MD;
global   Q1 Q2 Q3 Q4 U1 U2 U3 U4;
global   Q1p Q2p Q3p Q4p U1p U2p U3p U4p;
global   DEGtoRAD RADtoDEG COEF RHS SolutionToAxEqualsB;
global   TINITIAL TFINAL INTEGSTP PRINTINT ABSERR RELERR;

%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
LA                              =  2;                      % UNITS               Constant
LB                              =  2;                      % UNITS               Constant
LC                              =  2;                      % UNITS               Constant
LD                              =  2;                      % UNITS               Constant
MB                              =  5;                      % UNITS               Constant
MC                              =  5;                      % UNITS               Constant
MD                              =  5;                      % UNITS               Constant

Q1                              =  0;                      % UNITS               Initial Value
Q2                              =  0;                      % UNITS               Initial Value
Q3                              =  0;                      % UNITS               Initial Value
Q4                              =  0;                      % UNITS               Initial Value
U1                              =  1;                      % UNITS               Initial Value
U2                              =  1;                      % UNITS               Initial Value
U3                              =  1;                      % UNITS               Initial Value
U4                              =  1;                      % UNITS               Initial Value

TINITIAL                        =  0.0;                    % UNITS               Initial Time
TFINAL                          =  1.0;                    % UNITS               Final Time
INTEGSTP                        =  0.1;                    % UNITS               Integration Step
PRINTINT                        =  1;                      % Positive Integer    Print-Integer
ABSERR                          =  1.0E-08;                %                     Absolute Error
RELERR                          =  1.0E-07 ;               %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
Pi       = 3.141592653589793;
DEGtoRAD = Pi/180.0;
RADtoDEG = 180.0/Pi;

% Reserve space and initialize matrices
COEF = zeros(4,4);
RHS = zeros(1,4);

% Evaluate constants
% Set the initial values of the states
VAR(1) = Q1;
VAR(2) = Q2;
VAR(3) = Q3;
VAR(4) = Q4;
VAR(5) = U1;
VAR(6) = U2;
VAR(7) = U3;
VAR(8) = U4;



%===========================================================================
function OpenOutputFilesAndWriteHeadings
FileIdentifier = fopen('robot981.1', 'wt');   if( FileIdentifier == -1 ) error('Error: unable to open file robot981.1'); end
fprintf( 1,             '%%       T             Q1             Q2             Q3             Q4       DOT(P_NO_AO>,N DOT(P_NO_AO>,N DOT(P_NO_AO>,N DOT(P_AO_BO>,N DOT(P_AO_BO>,N DOT(P_AO_BO>,N DOT(P_BO_CO>,N DOT(P_BO_CO>,N DOT(P_BO_CO>,N DOT(P_CO_DO>,N DOT(P_CO_DO>,N DOT(P_CO_DO>,N\n' );
fprintf( 1,             '%%    (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
fprintf(FileIdentifier, '%% FILE: robot981.1\n%%\n' );
fprintf(FileIdentifier, '%%       T             Q1             Q2             Q3             Q4       DOT(P_NO_AO>,N DOT(P_NO_AO>,N DOT(P_NO_AO>,N DOT(P_AO_BO>,N DOT(P_AO_BO>,N DOT(P_AO_BO>,N DOT(P_BO_CO>,N DOT(P_BO_CO>,N DOT(P_BO_CO>,N DOT(P_CO_DO>,N DOT(P_CO_DO>,N DOT(P_CO_DO>,N\n' );
fprintf(FileIdentifier, '%%    (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );



%===========================================================================
% Main driver loop for numerical integration of differential equations
%===========================================================================
function SolveOrdinaryDifferentialEquations
global   LA LB LC LD MB MC MD;
global   Q1 Q2 Q3 Q4 U1 U2 U3 U4;
global   Q1p Q2p Q3p Q4p U1p U2p U3p U4p;
global   DEGtoRAD RADtoDEG COEF RHS SolutionToAxEqualsB;
global   TINITIAL TFINAL INTEGSTP PRINTINT ABSERR RELERR;

OpenOutputFilesAndWriteHeadings
VAR = ReadUserInput;
OdeMatlabOptions = odeset( 'RelTol',RELERR, 'AbsTol',ABSERR, 'MaxStep',INTEGSTP );
T = TINITIAL;
PrintCounter = 0;
mdlDerivatives(T,VAR,0);
while 1,
  if( TFINAL>=TINITIAL & T+0.01*INTEGSTP>=TFINAL ) PrintCounter = -1; end
  if( TFINAL<=TINITIAL & T+0.01*INTEGSTP<=TFINAL ) PrintCounter = -1; end
  if( PrintCounter <= 0.01 ),
     mdlOutputs(T,VAR,0);
     if( PrintCounter == -1 ) break; end
     PrintCounter = PRINTINT;
  end
  [TimeOdeArray,VarOdeArray] = ode45( @mdlDerivatives, [T T+INTEGSTP], VAR, OdeMatlabOptions, 0 );
  TimeAtEndOfArray = TimeOdeArray( length(TimeOdeArray) );
  if( abs(TimeAtEndOfArray - (T+INTEGSTP) ) >= abs(0.001*INTEGSTP) ) warning('numerical integration failed'); break;  end
  T = TimeAtEndOfArray;
  VAR = VarOdeArray( length(TimeOdeArray), : );
  PrintCounter = PrintCounter - 1;
end
mdlTerminate(T,VAR,0);



%===========================================================================
% mdlDerivatives: Calculates and returns the derivatives of the continuous states
%===========================================================================
function sys = mdlDerivatives(T,VAR,u)
global   LA LB LC LD MB MC MD;
global   Q1 Q2 Q3 Q4 U1 U2 U3 U4;
global   Q1p Q2p Q3p Q4p U1p U2p U3p U4p;
global   DEGtoRAD RADtoDEG COEF RHS SolutionToAxEqualsB;
global   TINITIAL TFINAL INTEGSTP PRINTINT ABSERR RELERR;

% Update variables after integration step
Q1 = VAR(1);
Q2 = VAR(2);
Q3 = VAR(3);
Q4 = VAR(4);
U1 = VAR(5);
U2 = VAR(6);
U3 = VAR(7);
U4 = VAR(8);
Q1p = U1;
Q2p = U2;
Q3p = U3;
Q4p = U4;

COEF(1,1) = -0.3333333333333333 - 0.25*MB*LB^2*sin(Q2)^2 - 0.25*MC*(4*LB^2*sin(Q2)^2+LC^2*cos(Q2+Q3)^2+4*LB*LC*sin(Q2)*cos(Q2+Q3)) - 0.25*MD*(4*LB^2*sin(Q2)^2+LD^2*cos(Q2+Q3)^2+4*LC*LD*cos(Q2+Q3)^2+4*LC^2*cos(Q2+Q3)^2+4*LB*LD*sin(Q2)*cos(Q2+Q3)+8*LB*LC*sin(Q2)*cos(Q2+Q3));
COEF(1,2) = 0;
COEF(1,3) = 0;
COEF(1,4) = 0.08333333333333333*sin(Q2+Q3);
COEF(2,1) = 0;
COEF(2,2) = -0.25 - 0.25*MB*LB^2 - 0.25*MC*(LC^2+4*LB^2-4*LB*LC*sin(Q3)) - 0.25*MD*(LD^2+4*LB^2+4*LC*LD+4*LC^2-8*LB*LC*sin(Q3)-4*LB*LD*sin(Q3));
COEF(2,3) = -0.1666666666666667 - 0.25*LC*MC*(LC-2*LB*sin(Q3)) - 0.25*MD*(LD+2*LC)*(LD+2*LC-2*LB*sin(Q3));
COEF(2,4) = 0;
COEF(3,1) = 0;
COEF(3,2) = -0.1666666666666667 - 0.25*LC*MC*(LC-2*LB*sin(Q3)) - 0.25*MD*(LD+2*LC)*(LD+2*LC-2*LB*sin(Q3));
COEF(3,3) = -0.1666666666666667 - 0.25*MC*LC^2 - 0.25*MD*(LD+2*LC)^2;
COEF(3,4) = 0;
COEF(4,1) = 0.08333333333333333*sin(Q2+Q3);
COEF(4,2) = 0;
COEF(4,3) = 0;
COEF(4,4) = -0.08333333333333333;
RHS(1) = -1 + sin(Q2+Q3) + 0.5*MB*LB^2*sin(Q2)*cos(Q2)*U1*U2 + 0.08333333333333333*sin(Q2+Q3)*cos(Q2+Q3)*U1*(U2+U3) + 0.08333333333333333*sin(Q4)*cos(Q2+Q3)*(cos(Q4)*cos(Q2+Q3)*U1*U4-sin(Q4)*U2*U4-sin(Q4)*U3*U4-sin(Q4)*sin(Q2+Q3)*U1*(U2+U3)) + 0.5*MC*U1*(4*LB^2*sin(Q2)*cos(Q2)*U2+2*LB*LC*cos(Q2)*cos(Q2+Q3)*U2-2*LB*LC*sin(Q2)*sin(Q2+Q3)*(U2+U3)-LC^2*sin(Q2+Q3)*cos(Q2+Q3)*(U2+U3)) + 0.25*MD*U1*(8*LB^2*sin(Q2)*cos(Q2)*U2+4*LB*LD*cos(Q2)*cos(Q2+Q3)*U2+8*LB*LC*cos(Q2)*cos(Q2+Q3)*U2-2*LB*LD*sin(Q2)*sin(Q2+Q3)*U2-2*LC*LD*sin(Q2+Q3)*cos(Q2+Q3)*U2-LD^2*sin(Q2+Q3)*cos(Q2+Q3)*U2-2*LB*LD*sin(Q2)*sin(Q2+Q3)*(U2+U3)-2*LC*LD*sin(Q2+Q3)*cos(Q2+Q3)*(U2+U3)-LD^2*sin(Q2+Q3)*cos(Q2+Q3)*(U2+U3)-2*LB*sin(Q2)*sin(Q2+Q3)*(LD*U3+4*LC*(U2+U3))-2*LC*sin(Q2+Q3)*cos(Q2+Q3)*(LD*U3+4*LC*(U2+U3))-LD*sin(Q2+Q3)*cos(Q2+Q3)*(LD*U3+4*LC*(U2+U3))) - 0.08333333333333333*cos(Q4)*cos(Q2+Q3)*(cos(Q4)*U2*U4+cos(Q4)*U3*U4+sin(Q4)*cos(Q2+Q3)*U1*U4+cos(Q4)*sin(Q2+Q3)*U1*(U2+U3));
RHS(2) = 0.08333333333333333*sin(Q4)*(cos(Q4)*U2*U4+cos(Q4)*U3*U4+sin(Q4)*cos(Q2+Q3)*U1*U4+cos(Q4)*sin(Q2+Q3)*U1*(U2+U3)) + 0.08333333333333333*cos(Q4)*(cos(Q4)*cos(Q2+Q3)*U1*U4-sin(Q4)*U2*U4-sin(Q4)*U3*U4-sin(Q4)*sin(Q2+Q3)*U1*(U2+U3)) + 0.25*MD*(2*LB*LD*cos(Q3)*U2^2+4*LB*LC*cos(Q3)*U2^2+2*LB*LD*sin(Q3)*sin(Q4)*cos(Q4)*sin(Q2+Q3)*U1*U2+LD*sin(Q4)*cos(Q4)*sin(Q2+Q3)*(LD+2*LC-2*LB*sin(Q3))*U1*U2+LD*sin(Q2+Q3)*(2*LB*sin(Q2)+LD*cos(Q2+Q3)+2*LC*cos(Q2+Q3))*U1^2+2*LC*sin(Q2+Q3)*(2*LB*sin(Q2)+LD*cos(Q2+Q3)+2*LC*cos(Q2+Q3))*U1^2-2*LC*LD*sin(Q4)*cos(Q4)*sin(Q2+Q3)*U1*U2-LD^2*sin(Q4)*cos(Q4)*sin(Q2+Q3)*U1*U2-2*LB*cos(Q3)*(U2+U3)*(LD*U3+2*LC*(U2+U3))-2*LB*cos(Q2)*(2*LB*sin(Q2)+LD*cos(Q2+Q3)+2*LC*cos(Q2+Q3))*U1^2-2*LB*LD*cos(Q3)*U2*(cos(Q4)*(cos(Q4)*U2+cos(Q4)*U3+sin(Q4)*cos(Q2+Q3)*U1)+sin(Q4)*(sin(Q4)*U2+sin(Q4)*U3-cos(Q4)*cos(Q2+Q3)*U1))) - 0.25*MB*LB^2*sin(Q2)*cos(Q2)*U1^2 - 0.25*MC*(2*LB*LC*cos(Q3)*(U2+U3)^2+2*LB*cos(Q2)*(2*LB*sin(Q2)+LC*cos(Q2+Q3))*U1^2-2*LB*LC*cos(Q3)*U2^2-LC*sin(Q2+Q3)*(2*LB*sin(Q2)+LC*cos(Q2+Q3))*U1^2);
RHS(3) = 0.25*LC*MC*(2*LB*cos(Q3)*U2^2+sin(Q2+Q3)*(2*LB*sin(Q2)+LC*cos(Q2+Q3))*U1^2) + 0.08333333333333333*sin(Q4)*(cos(Q4)*U2*U4+cos(Q4)*U3*U4+sin(Q4)*cos(Q2+Q3)*U1*U4+cos(Q4)*sin(Q2+Q3)*U1*(U2+U3)) + 0.08333333333333333*cos(Q4)*(cos(Q4)*cos(Q2+Q3)*U1*U4-sin(Q4)*U2*U4-sin(Q4)*U3*U4-sin(Q4)*sin(Q2+Q3)*U1*(U2+U3)) + 0.25*MD*(LD+2*LC)*(2*LB*cos(Q3)*U2^2+sin(Q2+Q3)*(2*LB*sin(Q2)+LD*cos(Q2+Q3)+2*LC*cos(Q2+Q3))*U1^2);
RHS(4) = -1 - 0.08333333333333333*cos(Q2+Q3)*U1*(U2+U3);
SolutionToAxEqualsB = COEF\RHS';

% Update variables after uncoupling equations
U1p = SolutionToAxEqualsB(1);
U2p = SolutionToAxEqualsB(2);
U3p = SolutionToAxEqualsB(3);
U4p = SolutionToAxEqualsB(4);

% Update derivative array prior to integration step
VARp(1) = Q1p;
VARp(2) = Q2p;
VARp(3) = Q3p;
VARp(4) = Q4p;
VARp(5) = U1p;
VARp(6) = U2p;
VARp(7) = U3p;
VARp(8) = U4p;

sys = VARp';



%===========================================================================
% mdlOutputs: Calculates and return the outputs
%===========================================================================
function Output = mdlOutputs(T,VAR,u)
global   LA LB LC LD MB MC MD;
global   Q1 Q2 Q3 Q4 U1 U2 U3 U4;
global   Q1p Q2p Q3p Q4p U1p U2p U3p U4p;
global   DEGtoRAD RADtoDEG COEF RHS SolutionToAxEqualsB;
global   TINITIAL TFINAL INTEGSTP PRINTINT ABSERR RELERR;

% Evaluate output quantities
Output(1)=T;  Output(2)=Q1;  Output(3)=Q2;  Output(4)=Q3;  Output(5)=Q4;  Output(6)=0.0;  Output(7)=0.0;  Output(8)=(0.5*LA);  Output(9)=(0.5*LB*sin(Q2)*cos(Q1));  Output(10)=(0.5*LB*sin(Q1)*sin(Q2));  Output(11)=(0.5*LA+0.5*LB*cos(Q2));  Output(12)=(0.5*LB*sin(Q2)*cos(Q1)+0.5*LC*cos(Q1)*cos(Q2+Q3));  Output(13)=(0.5*LB*sin(Q1)*sin(Q2)+0.5*LC*sin(Q1)*cos(Q2+Q3));  Output(14)=(0.5*LB*cos(Q2)-0.5*LC*sin(Q2+Q3));  Output(15)=(0.5*LC*cos(Q1)*cos(Q2+Q3)+0.5*LD*cos(Q1)*cos(Q2+Q3));  Output(16)=(0.5*LC*sin(Q1)*cos(Q2+Q3)+0.5*LD*sin(Q1)*cos(Q2+Q3));  Output(17)=(-0.5*LC*sin(Q2+Q3)-0.5*LD*sin(Q2+Q3));
FileIdentifier = fopen('all');
WriteOutput( 1,                 Output(1:17) );
WriteOutput( FileIdentifier(1), Output(1:17) );



%===========================================================================
function WriteOutput( fileIdentifier, Output )
numberOfOutputQuantities = length( Output );
if numberOfOutputQuantities > 0,
  for i=1:numberOfOutputQuantities,
    fprintf( fileIdentifier, ' %- 14.6E', Output(i) );
  end
  fprintf( fileIdentifier, '\n' );
end



%===========================================================================
% mdlTerminate: Perform end of simulation tasks and set sys=[]
%===========================================================================
function sys = mdlTerminate(T,VAR,u)
FileIdentifier = fopen('all');
fclose( FileIdentifier(1) );
fprintf( 1, '\n Output is in the file robot981.1\n' );
fprintf( 1, '\n To load and plot columns 1 and 2 with a solid line and columns 1 and 3 with a dashed line, enter:\n' );
fprintf( 1, '    someName = load( ''robot981.1'' );\n' );
fprintf( 1, '    plot( someName(:,1), someName(:,2), ''-'', someName(:,1), someName(:,3), ''--'' )\n\n' );
sys = [];



%===========================================================================
% Sfunction: System/Simulink function from standard template
%===========================================================================
function [sys,x0,str,ts] = Sfunction(t,x,u,flag)
switch flag,
  case 0,  [sys,x0,str,ts] = mdlInitializeSizes;    % Initialization of sys, initial state x0, state ordering string str, and sample times ts
  case 1,  sys = mdlDerivatives(t,x,u);             % Calculate the derivatives of continuous states and store them in sys
  case 2,  sys = mdlUpdate(t,x,u);                  % Update discrete states x(n+1) in sys
  case 3,  sys = mdlOutputs(t,x,u);                 % Calculate outputs in sys
  case 4,  sys = mdlGetTimeOfNextVarHit(t,x,u);     % Return next sample time for variable-step in sys
  case 9,  sys = mdlTerminate(t,x,u);               % Perform end of simulation tasks and set sys=[]
  otherwise error(['Unhandled flag = ',num2str(flag)]);
end



%===========================================================================
% mdlInitializeSizes: Return the sizes, initial state VAR, and sample times ts
%===========================================================================
function [sys,VAR,stateOrderingStrings,timeSampling] = mdlInitializeSizes
sizes = simsizes;             % Call simsizes to create a sizes structure
sizes.NumContStates  = 8;     % sys(1) is the number of continuous states
sizes.NumDiscStates  = 0;     % sys(2) is the number of discrete states
sizes.NumOutputs     = 17;    % sys(3) is the number of outputs
sizes.NumInputs      = 0;     % sys(4) is the number of inputs
sizes.DirFeedthrough = 1;     % sys(6) is 1, and allows for the output to be a function of the input
sizes.NumSampleTimes = 1;     % sys(7) is the number of samples times (the number of rows in ts)
sys = simsizes(sizes);        % Convert it to a sizes array
stateOrderingStrings = [];
timeSampling         = [0 0]; % m-by-2 matrix containing the sample times
OpenOutputFilesAndWriteHeadings
VAR = ReadUserInput
