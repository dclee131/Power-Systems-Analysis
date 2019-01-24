% written by Dongchan Lee
% Date: November, 2015
clear; close all;

%% System Data
%br_no fr to  R(PU)    X(pu)
LD=[1  1  2   0.04997  0.06644;
    2  2  3   0.02332  0.03310 ;
    3  1  4   0.04664  0.06201;
    4  1  5   0.02082  0.02768;
    5  5  6   0.02500  0.03322;
    6  1  7   0.02665  0.03543;
    7  7  8   0.02748  0.03654;
    8  1  9   0.03331  0.04430;
    9  1  10  0.02082  0.02768;
    10 2  11  0.02082  0.02768];

% bus no     activepower   reactivepower
BD=[1          0              0
    2         1.22          0.916
    3         0.032         0.024
    4         0.778         0.584
    5         0.673         0.595
    6         1.22          0.916
    7         0.0488        0.0366
    8         0.956         0.717
    9         0.698         0.523
    10        1.265         0.949
    11        0.265         0.0949];

num_bus=size(BD,1);
num_line=size(LD,1);
S=complex(BD(:,2),BD(:,3)); % Complex power injection
Z=complex(LD(:,4),LD(:,5));% branch impedance
line_from=LD(:,2);
line_to=LD(:,3);

A=zeros(num_bus,num_line); % A*I=IB where I is current injection and IB is branch current
for i=1:num_line
    A(line_from(i),i)=1;
    A(line_to(i),i)=-1;
end

%% Simulation parameters
V=ones(num_bus,1);% Initialize voltage
tol=1e-3; % Set the tolerence
itr_max=100; % Maximum iteration

%% Forward/Backward Method
% Refer to "A compensation-based power flow method for weakly meshed
% distribution and transmission networks" by D.Shirmohammadi et. al.
A_nonslack=A;
A_nonslack(1,:)=[];

i=1; S_mismatch=tol+1;
while ( S_mismatch(end) > tol && i < itr_max )
    I=conj(S./V); % Compute the injection current
    % backward sweep
    B=[zeros(num_line,1) inv(A_nonslack)];
    IB=B*I;  % caculate branch current
    
    % forward sweep
    V=ones(num_bus,1)-B'*diag(Z)*IB; % update voltage
    S_mismatch(i)=max(abs(S-V.*conj(I))); % Calculate S mismatch for terminal criteria
    V_hist(:,i)=V;
    i=i+1;
end

figure; % Plot convergence
plot(S_mismatch)
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('iteration'); ylabel('S mismatch')

%% Obtain every path from root to the leaves for plot
temp=1:num_line;
root=BD(1,1);
leaf=find(sum(abs(A),2)==1);
br_id=zeros(size(leaf,1),num_bus); % branch identification
br_id(:,1)=leaf;
temp(find(sum(A(leaf,temp))~=0))=[];

while ~isempty(temp)
    leaf_temp=find(sum(abs(A(:,temp)),2)==1);
    leaf_temp(find(leaf_temp==root))=[]; % if it's root delete
    
    for i=leaf_temp'
        cnt_bus=find(sum(abs(A(:,find(A(i,:)~=0))),2)~=0); % find lines that are connected to i
        cnt_bus(find(cnt_bus==root))=[]; % delete root
        cnt_bus(find(cnt_bus==i))=[]; % delete itself
        for j=cnt_bus'
            for k=find(br_id(:,1)==j)
                br_id(k,:)=[i br_id(k,1:end-1)];
            end
        end
    end
    temp(find(sum(A(leaf_temp,temp))~=0))=[];
    if isempty(temp) && sum(br_id(:,1)~=root)~=0
        br_id=[root*ones(size(leaf,1),1) br_id(:,1:end-1)];
    end
end
br_id(:,~any(br_id,1))=[]; % Delete auxilary zero columns

figure;
subplot(1,2,1) % voltage magnitude along every path
grid on;
hold all;
for i=1:size(br_id,1)
    bus_path=br_id(i,:);
    bus_path(~any(bus_path,1))=[];
    plot(1:size(bus_path,2),abs(V(bus_path)),'.-')
end

subplot(1,2,2) % voltage vectors along every path
grid on;
for i=1:size(br_id,1)
    bus_path=br_id(i,:);
    bus_path(~any(bus_path,1))=[];
    polar(angle(V(bus_path)),abs(V(bus_path)),'o:')
    hold on;
end

