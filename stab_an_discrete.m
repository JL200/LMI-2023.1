%%% Class 1 - TASK
%%% João Luiz
clear all
close all
clc

%% discrete system:

% OL:
%A'*P*A-P<0
%P>0

Ts = 0.1;

num = 1;
% den = [1 -0.1 1];            % unstable
den = [1 0.1 1];            % estable

G = tf(num,den);
Gd = c2d(G,Ts,'zoh')
% step(Gd)

[b,a] = tfdata(Gd,'v')

[A,B,C,D] = tf2ss(b,a)
sys = ss(A,B,C,D,Ts);

figure;
step(G)
hold on
grid on
step(sys,'r')

%% CL system:
% Q*A'+N'*B'<0
% Q>0
[x,y] = size(B);

Q = sdpvar(length(A));
N = sdpvar(y,x);

lmi1 = [[-Q Q*A'+N'*B'; A*Q+B*N -Q]<=0];
lmi2 = [Q>=0];

LMI = [lmi1,lmi2];

optimize(LMI)

Qs = value(Q);
Ns = value(N);

eig(A)
eig(Q*A'+N'*B')
eig(Qs)

K=Ns*inv(Qs);

checkset(LMI)

figure;
sys = ss(A+B*K,B,C,D,Ts)
step(sys)
grid on