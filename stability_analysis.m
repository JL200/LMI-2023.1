%%% AULA 1
%%% João Luiz
clear all
close all
clc

%% continuous system:
%A'P + PA < 0
%P > 0

num = 1;
den = [1 -0.1 1];            % unstable
% den = [1 0.1 1];            % estable

[A,B,C,D] = tf2ss(num,den);
sistem = ss(A,B,C,D);

step(sistem)

%% stability analysis

P = sdpvar(2,2)
sdisplay(P)

lmi1 = [A'*P+P*A<=0]
lmi2 = [P>=0]

LMI = [lmi1,lmi2];

optimize(LMI);

Ps = value(P);
eig(A)

checkset(LMI)

eig(A'*Ps+Ps*A)

%% discrete system:


