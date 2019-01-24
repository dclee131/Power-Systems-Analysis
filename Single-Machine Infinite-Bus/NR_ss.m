function x = NR_ss(Y,S_inj,idx_pq,slack_bus)
% Newton Raphson to solve steady state
max_iter=25;
num_bus=size(Y,1);
v_mag=1+0.01*rand(num_bus,1); v_mag(slack_bus)=1;
theta=0.01*rand(num_bus,1); theta(slack_bus)=0;
idx_nslack=[1:slack_bus-1 slack_bus+1:num_bus];
for iter=1:max_iter
    v_cpx=v_mag.*cos(theta)+1i*v_mag.*sin(theta); % complex voltage
    S_bal=v_cpx.*conj(Y*v_cpx)-S_inj; % Apparent power imbalance
    f=[real(S_bal(idx_nslack)); imag(S_bal(idx_pq))];
    J1=real(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(1i*v_cpx).'))+diag(1i*v_cpx)*diag(conj(Y*v_cpx)));
    J2=real(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(v_cpx./v_mag).'))+diag(v_cpx./v_mag)*diag(conj(Y*v_cpx)));
    J3=imag(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(1i*v_cpx).'))+diag(1i*v_cpx)*diag(conj(Y*v_cpx)));
    J4=imag(diag(v_cpx)*conj(Y.*(ones(num_bus,1)*(v_cpx./v_mag).'))+diag(v_cpx./v_mag)*diag(conj(Y*v_cpx)));
    J =[J1(idx_nslack,idx_nslack) J2(idx_nslack,idx_pq); J3(idx_pq,idx_nslack) J4(idx_pq,idx_pq)];
    dx=-J\f;
    
    nf(iter)=norm(f); ndx(iter)=norm(dx);
    theta(idx_nslack)=theta(idx_nslack)+dx(1:num_bus-1);
    v_mag(idx_pq)=v_mag(idx_pq)+dx(num_bus:end);
    if (nf(iter)<1e-8)&&(ndx(iter)<1e-8); break; end;
end
x=[theta; v_mag];
if iter==max_iter; disp('Newton Raphson did not converge!'); x=inf; end
%hold all; subplot(2,1,1); plot(dx); subplot(2,1,2); plot(nf)
