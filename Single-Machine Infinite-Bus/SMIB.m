% Written by Dongchan Lee, August 2017
clear; close all;

%% Select options
t_end=20;
del_t=0.05;

%% Set up variables
slack_bus=1;
num_bus=2;
num_line=1;
num_gen=1;
num_load=1;

eye_bus=eye(num_bus);
idx_delta=1:num_gen; idx_omega=(1:num_gen)+num_gen; 
idx_eq=(1:num_gen)+2*num_gen; idx_ef=(1:num_gen)+3*num_gen;
idx_id=(1:num_gen)+4*num_gen; idx_iq=(1:num_gen)+5*num_gen;
idx_Vd=(1)+6*num_gen; idx_Vq=(1)+6*num_gen+1; 
idx_dvar=1:4*num_gen; idx_avar=4*num_gen+1:6*num_gen+2*num_bus;

idx_gen=2;
idx_load=[];

M_gen=8;
D_gen=2;
ra=0;
xd=0.5;
xd_p=0.05;
xq=0.4;
xq_p=0.1;
Td=6;
TA=0.2;
KA=30;
line_frto=[1 2];
Z_line=0.03+1i*0.1;
E=eye_bus(line_frto(:,1),:)-eye_bus(line_frto(:,2),:);
Y_pref=E'*diag(Z_line.^-1)*E;

Pgen=zeros(num_bus,1); Sload=zeros(num_bus,1);
Pgen(2)=1;
Sload(1)=(1+1i*0.2);
S_inj=(Pgen-Sload);

Y=Y_pref;

%% Initialize data
% Initialize after conversion to impedence load
S_inj(idx_load)=0;
theta=NR_ss(Y,S_inj,idx_load,slack_bus); % Power flow
V_eq=theta(num_bus+1:end).*(cos(theta(1:num_bus))+1i*sin(theta(1:num_bus)));
I_eq=Y*V_eq;
omega_eq=zeros(num_gen,1);
delta_eq=angle(V_eq(idx_gen)+(ra+1i*xq_p).*I_eq(idx_gen));
v_eq=V_eq(idx_gen).*exp(-1i*(delta_eq-pi/2)); vd_eq=real(v_eq); vq_eq=imag(v_eq);
i_eq=I_eq(idx_gen).*exp(-1i*(delta_eq-pi/2)); id_eq=real(i_eq); iq_eq=imag(i_eq);
eq_eq=vq_eq+xd_p.*id_eq+ra.*iq_eq;
ef_eq=eq_eq+(xd-xd_p).*id_eq;
Pm=(vd_eq+ra.*id_eq).*id_eq+(vq_eq+ra.*iq_eq).*iq_eq;
x_eq=[delta_eq; omega_eq; eq_eq; ef_eq; id_eq; iq_eq; real(V_eq(2)); imag(V_eq(2))];
delta=delta_eq; omega=omega_eq; eq=eq_eq; id=id_eq; iq=iq_eq; vd_grid=real(V_eq(2)); vq_grid=imag(V_eq(2)); % for debugging
v_ref=abs(V_eq(2))^2+ef_eq/KA;

%% Define function
ff=@(delta,omega,eq,ef,id,iq,vd_grid,vq_grid) [omega;
       (-D_gen.*omega+Pm-eq.*iq-(xq_p-xd_p).*id.*iq)./M_gen;
       (ef-eq-(xd-xd_p).*id)./Td;
       (-ef+KA*(v_ref-vd_grid^2-vq_grid^2))/TA];
f_x=@(x) ff(x(idx_delta),x(idx_omega),x(idx_eq),x(idx_ef),x(idx_id),x(idx_iq),x(idx_Vd),x(idx_Vq));

JJ_fx=@(delta,omega,eq,ef,id,iq,vd_grid,vq_grid) [zeros(num_gen) eye(num_gen) zeros(num_gen,2*num_gen);
    zeros(num_gen) diag(-D_gen./M_gen) diag(-iq./M_gen) zeros(num_gen);
    zeros(num_gen,2*num_gen) diag(-1./Td) diag(1./Td)
    zeros(num_gen,2*num_gen) diag(1./TA) zeros(num_gen)];

JJ_fy=@(delta,omega,eq,id,iq,vd_grid,vq_grid) [zeros(num_gen,2*num_gen+2);
    diag((-(xq_p-xd_p).*iq)./M_gen) diag((-eq-(xq_p-xd_p).*id)./M_gen) zeros(num_gen,2);
    diag(-(xd-xd_p)./Td) zeros(num_gen) zeros(num_gen,2)
    zeros(num_gen,2*num_gen) diag(-2*vd_grid./TA) diag(-2*vq_grid./TA)];

gg=@(delta,eq,id,iq,vd_grid,vq_grid) [-ra.*id+xq_p.*iq-vd_grid.*sin(delta)+vq_grid.*cos(delta);
    eq-ra.*iq-xd_p.*id-vd_grid.*cos(delta)-vq_grid.*sin(delta);
    real(Y(2,:))*[1; vd_grid]-imag(Y(2,:))*[0; vq_grid]-(id.*sin(delta)+iq.*cos(delta));
    imag(Y(2,:))*[1; vd_grid]+real(Y(2,:))*[0; vq_grid]-(-id.*cos(delta)+iq.*sin(delta))];

g_x=@(x) gg(x(idx_delta),x(idx_eq),x(idx_id),x(idx_iq),x(idx_Vd),x(idx_Vq));

JJ_gx=@(delta,omega,eq,id,iq,vd_grid,vq_grid) [diag(-vd_grid.*cos(delta)-vq_grid.*sin(delta)) zeros(num_gen,3*num_gen);
    diag(vd_grid.*sin(delta)-vq_grid.*cos(delta)) zeros(num_gen,num_gen) eye(num_gen) zeros(num_gen);
    -diag(id.*cos(delta)-iq.*sin(delta)) zeros(1,3*num_gen);
    -diag(id.*sin(delta)+iq.*cos(delta)) zeros(1,3*num_gen)];

JJ_gy=@(delta,omega,eq,id,iq,vd_grid,vq_grid) [diag(-ra) diag(xq_p) diag(-sin(delta)) diag(cos (delta));
    diag(-xd_p) diag(-ra) diag(-cos(delta)) diag(-sin(delta));
    -diag(sin(delta)) -diag(cos(delta)) real(Y(2,2)) -imag(Y(2,2));
    diag(cos(delta)) -diag(sin(delta)) imag(Y(2,2)) real(Y(2,2))];

J_fx=@(x) JJ_fx(x(idx_delta),x(idx_omega),x(idx_eq),x(idx_id),x(idx_iq),x(idx_Vd),x(idx_Vq));
J_fy=@(x) JJ_fy(x(idx_delta),x(idx_omega),x(idx_eq),x(idx_id),x(idx_iq),x(idx_Vd),x(idx_Vq));
J_gx=@(x) JJ_gx(x(idx_delta),x(idx_omega),x(idx_eq),x(idx_id),x(idx_iq),x(idx_Vd),x(idx_Vq));
J_gy=@(x) JJ_gy(x(idx_delta),x(idx_omega),x(idx_eq),x(idx_id),x(idx_iq),x(idx_Vd),x(idx_Vq));
J_x=@(x) [J_fx(x) J_fy(x); J_gx(x) J_gy(x)];
fg_x=@(x) [f_x(x); g_x(x)];
J_gy=@(x) JJ_gy(x(idx_delta),x(idx_Vd),x(idx_Vq));

% Code for checking function and its Jacobian check by perturbation
%sum(fg_x(x_eq)>1e-10) % number of violations in f
% x=rand(6*num_gen+2*num_bus,1); J_check=[];
% for i=1:size(x,1)
%     x_unit=zeros(size(x,1),1); x_unit(i)=1;
%     J_check(:,i)=(fg_x(x+0.001*x_unit)-fg_x(x))/0.001;
% end
% J_fail=max(max(abs(J_x(x)-J_check)));

%% Simulation
t=1:del_t:t_end;
x_sim=x_eq;
x_sim(1)=x_sim(1)+0.1;

for t_idx=1:t_end/del_t
    f=@(x) [-x(idx_dvar)+x_sim(idx_dvar,t_idx)+del_t/2*(f_x(x_sim(:,t_idx))+f_x(x)); g_x(x)]; % Trapezoidal rule
    J=@(x) [-eye(size(idx_dvar,2))+del_t/2*J_fx(x) del_t/2*J_fy(x); J_gx(x) J_gy(x)];
    x_sim(:,t_idx+1)=NR(f,J,x_sim(:,t_idx));
    if sum(x_sim(:,t_idx+1)==inf); break; end
end

plot(x_sim(idx_dvar,:)')

%% Set up constants for linearized system at the equilibrium
r_line=real(Z_line); x_line=imag(Z_line);
Delta=r_line^2+(x_line+xq_p)*(x_line+xd_p);
k1=1/Delta*(iq_eq*(xq_p-xd_p)*(-r_line*cos(delta_eq)+(x_line+xq_p)*sin(delta_eq))+(eq_eq+(xq_p-xd_p)*id_eq)*((x_line+xd_p)*cos(delta_eq)+r_line*sin(delta_eq)));
k2=iq_eq+1/Delta*((xq_p-xd_p)*(x_line+xq_p)*iq_eq+r_line*eq_eq+r_line*(xq_p-xd_p)*id_eq);
k3=(1+(xd-xd_p)*(x_line+xq_p)/Delta)^-1;
k4=(xd-xd_p)/Delta*(-r_line*cos(delta_eq)+(x_line+xq_p)*sin(delta_eq));
vabs_eq=sqrt(vd_eq^2+vq_eq^2);
k5=1/Delta*(vd_eq/vabs_eq*xq_p*((x_line+xd_p)*cos(delta_eq)+r_line*sin(delta_eq))-vq_eq/vabs_eq*xd_p*(-r_line*cos(delta_eq)+(x_line+xq_p)*sin(delta_eq)));
k6=1/Delta*(vd_eq/vabs_eq*xq_p*r_line+vq_eq/vabs_eq*xd_p*(x_line+xq_p));

J_fx(x_eq)