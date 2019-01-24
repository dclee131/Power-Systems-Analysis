%% Written by Dongchan Lee, November 2016
%clear; close all;

%% Load case study and initialize
mpc=eval('case118'); % Imports data from MATPOWER
num_bus=size(mpc.bus,1);
num_line=size(mpc.branch,1);
line_frto=mpc.branch(:,1:2);
Z_line=mpc.branch(:,3)+1i*mpc.branch(:,4);
P_gen=zeros(num_bus,1); S_load=zeros(num_bus,1);
P_gen(mpc.gen(:,1))=mpc.gen(:,2)/mpc.baseMVA;
S_load(mpc.bus(:,1))=-(mpc.bus(:,3)+1i*mpc.bus(:,4))/mpc.baseMVA;
S_inj=P_gen+S_load;
idx_pq = find(mpc.bus(:,2)==1); % PQ Bus index

E=zeros(num_line,num_bus);
for i=1:num_line
    E(i,line_frto(i,1))=1;  E(i,line_frto(i,2))=-1;
end
Y=E'*diag(Z_line.^-1)*E; % Admittance matrix

%% Newton Raphson to solve steady state
max_iter=10;
v_mag=ones(num_bus,1);
theta(1)=0;
theta(2:num_bus,1)=imag(Y(2:num_bus,2:num_bus))\real(S_inj(2:num_bus));
for iter=1:max_iter
    v_cpx=v_mag.*cos(theta)+1i*v_mag.*sin(theta); % complex voltage
    S_bal=v_cpx.*conj(Y*v_cpx)-S_inj; % Apparent power imbalance
    f=[real(S_bal(2:end)); imag(S_bal(idx_pq))];
    J1=real(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(1i*v_cpx).'))+diag(1i*v_cpx)*diag(conj(Y*v_cpx)));
    J2=real(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(v_cpx./v_mag).'))+diag(v_cpx./v_mag)*diag(conj(Y*v_cpx)));
    J3=imag(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(1i*v_cpx).'))+diag(1i*v_cpx)*diag(conj(Y*v_cpx)));
    J4=imag(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(v_cpx./v_mag).'))+diag(v_cpx./v_mag)*diag(conj(Y*v_cpx)));
    J =[J1(2:end,2:end) J2(2:end,idx_pq); J3(idx_pq,2:end) J4(idx_pq,idx_pq)];
    dx=-J\f;
    
    nf(iter)=norm(f); ndx(iter)=norm(dx);
    theta(2:num_bus)=theta(2:num_bus)+dx(1:num_bus-1);
    v_mag(idx_pq)=v_mag(idx_pq)+dx(num_bus:end);
    if (nf(iter)<1e-8)&&(ndx(iter)<1e-8); break; end;
end
if iter==max_iter; disp('Newton Raphson did not converge!'); end

hold all;
plot(ndx,'b'); plot(nf,'r');
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('iteration'); ylabel('mismatch'); legend('|dx|','|f|')

