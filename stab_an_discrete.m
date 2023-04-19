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

[numd,dend] = tfdata(Gd,'v')

[A,B,C,D] = tf2ss(numd,dend)
sistem = ss(A,B,C,D);

step(sistem)

%% OL system:

P = sdpvar(length(A));
% sdisplay(P)

lmi1 = [A'*P*A-P<=0];
lmi2 = [P>=0];

LMI = [lmi1,lmi2];

optimize(LMI);

Ps = value(P);

eig(A)

checkset(LMI)

eig(A'*P*A-P)

checkset(LMI)

%% CL system:
% Q*A'+N'*B'<0
% Q>0
[x,y] = size(B);

Q = sdpvar(length(A));
N = sdpvar(y,x);

lmi1 = [Q*A'+N'*B'<=0];
lmi2 = [Q>=0];

LMI = [lmi1,lmi2];

optimize(LMI)

Qs = value(Q);
Ns = value(N);

eig(A)
eig(Q*A'+N'*B')
eig(Qs)

Q_aux = inv(Qs);
K=Ns*Q_aux;

checkset(LMI)

sys = ss(A-B*K,B,C-D*K,D)
step(sys)
grid on
