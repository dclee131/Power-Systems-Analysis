function x = NR(f_func,J_func,x0,accuracy,max_itr)
x=x0;
if nargin<4; accuracy=1e-6; end
if nargin<5; max_itr=30; end
for iter=1:max_itr % Newton loop
    J=J_func(x);
    f=f_func(x);
    dx=-J\f;
    nf(iter)=norm(f); ndx(iter)=norm(dx); % save norms for debugging
    if nf(iter) < accuracy && ndx(iter) < accuracy; break; end
    x=x+dx; % Update solution
    %x=max(x,-pi/2); x=min(x,pi/2);
end
if iter>=max_itr
    disp('NR did not converge!'); x=inf;
    hold all; subplot(2,1,1); plot(ndx); subplot(2,1,2); plot(nf)
end