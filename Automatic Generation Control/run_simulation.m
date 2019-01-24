% Written by Dongchan Lee
% Date: July 26, 2017

clear; close all;
H=5;% Inertia constant
S_B=10000;
w0=50; % Nominal frequency
J=2*H*S_B/w0^2; % Inertia
T_t=0; % Turbine parameter
D_L=1/200; W_0=0; % Load parameter
S=0.04; % Primary control parameter
C_p=0.17; T_N=120; % AGC parameter
B_freq=1/S; % Frequency bias

%%
[A,B,C,D]=linmod('one_area_model');
[b_cl,a_cl] = ss2tf(A,B(:,2),C,D(:,2));
H_cl=tf(b_cl,a_cl);

B_freq=0;
[A,B,C,D]=linmod('one_area_model');
[b,a] = ss2tf(A,B(:,2),C,D(:,2));
G=tf(b,a);
D_ctrl=tf([T_N*C_p 1],[T_N*C_p 0]);

figure; 
subplot(1,2,1); hold all; % step response
opt=stepDataOptions; opt.StepAmplitude=500; opt.InputOffset = 0;
t = (0:0.05:1800)';
w_AGC=step(H_cl,t, opt);
w_noAGC=step(G,t, opt);
plot(t,w_AGC+w0,'b-')
plot(t,w_noAGC+w0,'r--')
legend('with AGC','without AGC'); xlabel('time'); ylabel('f');

subplot(1,2,2);  % root locus
rlocus(G*D_ctrl)

%%
P_hat=1/S;
B_freq=1/S; % Frequency bias
[A,B,C,D]=linmod('two_area_model');
sys=ss(A,B(:,2),C,D(:,2));

figure; % step response
subplot(2,1,1); hold all;
opt=stepDataOptions; opt.StepAmplitude=10; opt.InputOffset = 0;
t = (0:0.05:2400)';
y_AGC=step(sys,t, opt);
plot(t,y_AGC(:,1:2)+w0)
legend('f1','f2'); xlabel('time'); ylabel('f');

subplot(2,1,2); hold all;
plot(t,y_AGC(:,3:5))
legend('\Delta P_{m1}','\Delta P_{m2}','\Delta P_{T12}'); xlabel('time'); ylabel('\Delta P');
