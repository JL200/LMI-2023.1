%%% João Luiz de Castro Pereira - UFC
%%% DC motor analysis
%%% LMI

clear all
close all
clc

%% DC motor especifications:
Leq = 0.1;
Req = 0.2;
J = 1
B1 = 0.1;
B2 = 1;
Ka = 1;

%% A:
A1 = [-Req/Leq 0 0;
    Ka/J -B1/J 0;
    0 0 1];
A2 = [-Req/Leq 0 0;
    Ka/J -B2/J 0;
    0 0 1];

% B:
B = [1/Leq; 0; 0];

% C:
C = [0 0 1];

% D:
D = 0;

%% 
[m,n] = size(B)
Q = sdpvar(n,n,'symmetric');
Y = sdpvar(n,m);

%% LMI

LMI1 = [Q>=0];
LMI2 = [(A1*Q+B*Y)+(A1*Q+B*Y)'<=0];
LMI3 = [(A2*Q+B*Y)+(A2*Q+B*Y)'<=0];
LMI = [LMI1,LMI2,LMI3];

optimize(LMI,-trace(Q));
checkset(LMI)

Y = value(Y)
Q = value(Q)

K = Y+inv(Q)

e1 = eig(A1+B*K)
e2 = eig(A2+B*K)