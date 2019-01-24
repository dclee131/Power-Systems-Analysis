%% Written by Dongchan Lee, July 2017
% We assume proportional control for AGC without considering interarea exchange error.
% So this is really just with primary frequency control.
clear; close all;

%% Load case study and initialize
mpc=eval('case39'); % Imports data from MATPOWER
num_bus=size(mpc.bus,1);
num_line=size(mpc.branch,1);
line_frto=mpc.branch(:,1:2);
Z_line=mpc.branch(:,3)+1i*mpc.branch(:,4);
P_gen=zeros(num_bus,1); S_load=zeros(num_bus,1);
P_gen(mpc.gen(:,1))=mpc.gen(:,2)/mpc.baseMVA;
S_load(mpc.bus(:,1))=-(mpc.bus(:,3)+1i*mpc.bus(:,4))/mpc.baseMVA;
S_inj=P_gen+S_load;
idx_pq = find(mpc.bus(:,2)==1); % PQ Bus index
S_prc=zeros(num_bus,1); k_AGC=zeros(num_bus,1);
S_prc(mpc.gen(:,1))=1/0.04;
k_AGC(mpc.gen(:,1))=1/0.0017;

E=zeros(num_line,num_bus);
for i=1:num_line
    E(i,line_frto(i,1))=1;  E(i,line_frto(i,2))=-1;
end
Y=E'*diag(Z_line.^-1)*E; % Admittance matrix

%% Newton Raphson to solve steady state
max_iter=100;
v_mag=ones(num_bus,1);
theta=[0; 0.01*rand(num_bus-1,1)];
omega=0;
for iter=1:max_iter
    v_cpx=v_mag.*cos(theta)+1i*v_mag.*sin(theta); % complex voltage
    S_bal=v_cpx.*conj(Y*v_cpx)-S_inj+(S_prc+k_AGC)*omega; % Apparent power imbalance
    f=[real(S_bal); imag(S_bal(idx_pq))];
    J11=real(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(1i*v_cpx).'))+diag(1i*v_cpx)*diag(conj(Y*v_cpx))); %dP/dtheta
    J12=real(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(v_cpx./v_mag).'))+diag(v_cpx./v_mag)*diag(conj(Y*v_cpx))); %dP/dV_mag
    J13=S_prc+k_AGC; %dP/domega
    J21=imag(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(1i*v_cpx).'))+diag(1i*v_cpx)*diag(conj(Y*v_cpx))); %dQ/dtheta
    J22=imag(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(v_cpx./v_mag).'))+diag(v_cpx./v_mag)*diag(conj(Y*v_cpx))); %dQ/dV_mag
    J23=zeros(num_bus,1);
    J =[J11 J12(:,idx_pq) J13; J21(idx_pq,:) J22(idx_pq,idx_pq) J23(idx_pq)];
    dx=-J\f;
    
    nf(iter)=norm(f); ndx(iter)=norm(dx);
    theta=theta+dx(1:num_bus);
    v_mag(idx_pq)=v_mag(idx_pq)+dx(num_bus+1:size(dx,1)-1);
    omega=omega+dx(end);
    if (nf(iter)<1e-8)&&(ndx(iter)<1e-8); break; end;
end
if iter==max_iter; disp('Newton Raphson did not converge!'); end

hold all;
plot(ndx,'b'); plot(nf,'r');
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('iteration'); ylabel('mismatch'); legend('|dx|','|f|')

