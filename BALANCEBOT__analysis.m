% (c) daniel ford, daniel.jb.ford@gmail.com

% controller analysis for Lego balancing robot

clear all;
close all;
close all hidden;

% constants
Tsamp = 0.001;  % sampling time
%Tsamp = 0.1;
Tend = 10;      % sim ending time, sec
mass = 0.8;     % kg (guesstimate)
L = 0.1651;     % m (6.5 inches) (measured)
g = 9.80665;    % m/s^2
b = 0.003;      % angular damping constant (WAG)
D2L = 5.5;      % conversion factor of degrees tilt to light sensor output
L2D = 1/D2L;    % inverse conversion factor   

% initial conditions
theta_IC = 3*pi/180;
thetaDot_IC = 0.0;
initial = [theta_IC thetaDot_IC]';

A = zeros(2);
B = A(:,1);
C = [1 0];
D = 0;

% 2x2 version
A(1,1:2) = [-b/(mass*L^2) g/L];
A(2,1) = 1;
B(1) = 1/(mass*L^2);

% generate system, decoupled version, discretized version
sys = ss(A,B,C,D);
sys_dec = My_Decouple(sys);
sys_z = c2d(sys_dec,Tsamp,'zoh');
[num,den] = ss2tf(sys_z.a,sys_z.b,sys_z.c,sys_z.d);

% choose poles

Pc = exp(Tsamp*[-0.0479;-0.0021]);
Po = exp(Tsamp*-.05);

% calc polynomials and Sylvester matrix
Dcl = conv([1 -Pc(1)],[1 -Pc(2)]);
Dob = [1 -Po]

s = conv(Dob,(Dcl-den));
Dstar = [s(4);s(3);s(2);s(1)];

sylv = [den(3)  0       num(3)    0;
        den(2)  den(3)  num(2)    num(3);
        1       den(2)  num(1)    num(2);
        0       1       0         num(1)];
        
C_fn = inv(sylv)*Dstar;

num0 = [C_fn(2) C_fn(1)]
den0 = Dob;
num1 = [C_fn(4) C_fn(3)]
den1 = Dob;

den0_num0 = den0 + num0;

% run simulation
LegoBalance_no_loop
set_param('LegoBalance_no_loop', 'SimulationCommand', 'start')
