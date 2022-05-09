clear all
close all
clc

%% values
A=[0 1; 880.87 0];
B=[0 -9.9453]';
C=[708.27 0];
D=[0];
x0_l=[1 1];
%if we want to compare the models you have to keep it fixed
%x0_f=1*rand(2, 6); 
x0_f=[0.1818 0.1455; 0.8693 0.5499; 0.8530 0.3510; ...
    0.2638 0.1361; 0.5797 0.1450; 0.6221 0.5132]';
%n=5; %noise
n=0; %no_noise
%% eigenvalues of A to change dynamics of reference signal
%K0=place(A, B, [-1 -2]);
%K0=place(A, B, [-1 0]);
%K0=acker(A, B, [0 0]);
K0=place(A, B, [-j +j]);
eig(A)
A=A-B*K0;
eig(A)

%% control variables
Adj=[0 0 0 0 0 0;
    2 0 0 0 0 0;
    0 3 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 3 0;];
DD=diag([0 2 3 1 1 3]);
% G=diag([1 0 0 0 0 0]);
G=diag([1 0 1 0 1 0]);
% Adj=zeros(6);
% DD=diag(zeros(1,6));
% G=diag(ones(1,6));
L=DD-Adj;
G=diag([1 0 0 0 0 0]);
G_ext=zeros(7, 7);
G_ext(2:end,2:end)=Adj;
G_ext(2:end,1)=diag(G);
lam_g=eig(L+G);

c=10*1/(2*min(real(lam_g)));
Q=diag([1 1]);
R=10;
P=are(A, B*inv(R)*B', Q);
K=R^-1*B'*P;
%% dist observer
Q=diag([1 1]);
R=1;
P=are(A', C'*inv(R)*C, Q);
F=P*C'*R^-1;

%% sim
open('homework_2_sim_dist_obs.slx');
out=sim('homework_2_sim_dist_obs.slx');
X=out.X.Data;
t=out.X.Time;
U=out.U.Data;

%Global disagreement error
delta1 = abs(X(:,3:2:length(X(1,:)))-X(:,1));
delta2 = abs(X(:,4:2:length(X(1,:)))-X(:,2));
figure(1)
hold on,grid on
plot(t, delta1, 'LineWidth', 1.5)
xlabel('time [s]')
ylabel('disagreement error [cm]')
title('Global disagreement error on x1')
legend('x1', 'x2', 'x3', 'x4', 'x5', 'x6')

figure(2)
hold on, grid on
plot(t, delta2, 'LineWidth', 1.5)
title('Global disagreement error on x2')
xlabel('time [s]')
ylabel('disagreement error [cm]')
legend('x1', 'x2', 'x3', 'x4', 'x5', 'x6')

figure(3)
hold on, grid on
plot(t, U, 'LineWidth', 1.5)
ylim([-4 4])
title('u')
xlabel('time [s]')
ylabel('command input [A]')
legend('u1', 'u2', 'u3', 'u4', 'u5', 'u6')

figure(4)
hold on, grid on
for i=1:2:14
    plot(t, X(:,i), 'LineWidth', 1.5)
end
title('x1')
xlabel('time [s]')
legend('x0', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6')

figure(5)
hold on, grid on
for i=2:2:14
    plot(t, X(:,i), 'LineWidth', 1.5)
end
title('x2')
xlabel('time [s]')
legend('x0', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6')
%%
% eig(L+G)
lam_g=eig(L+G);
for i = 1:length(lam_g)
    eig(A-c.*lam_g(i).*B*K)
end