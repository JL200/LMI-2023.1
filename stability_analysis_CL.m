%%% Class 1
%%% João Luiz
clear all
close all
clc

%% continuous system:
%QA'+AQ+Y'B'+BY<0
%Q>0

num = 1;
% den = [1 -0.1 1];            % unstable
den = [1 0.2 1];            % estable

[A,B,C,D] = tf2ss(num,den);
sistem = ss(A,B,C,D);



step(sistem)

%% CL system:

[x,y] = size(B)

Q = sdpvar(length(A));
Y = sdpvar(y,x);
% sdisplay(P)

lmi1 = [Q*A'+A*Q+Y'*B'+B*Y<=0];
lmi2 = [Q>=0];

LMI = [lmi1,lmi2];

optimize(LMI);

Qs = value(Q);
Ys = value(Y);

eig(A)

checkset(LMI)

eig(Q*A'+A*Q+Y'*B'+B*Y)

Q_aux = inv(Qs)
K = Ys*Q_aux;


checkset(LMI)