function [Theta_HC0,theta_bar0,eta_0,K,P,mu,H,rho_theta0,L_B,d_bar,c,c_max,L_cost,l_maxsteady,s_steady]=offlineComputation(Q,R,F,G,B_p,W,A_0,A_1,A_2,B_0,B_1,B_2,m,n,p,L_const,l_const)
options = optimset('Display','off');
%Compute inital parameter set Theta_HC0 (fixed)
theta_bar0= zeros(p,1);
eta_0=2;
Theta_HC0=theta_bar0+eta_0*B_p;
%Compute feedback K and terminal cost P
rho_PK = 0.75;%Hyperparameter
[P,K]=get_P_K(rho_PK,Q,R,A_0,A_1,A_2,B_0,B_1,B_2,Theta_HC0,m,n,p);
%Compute the offline computed polytope P
lambda=0.75;%Hyperparameter
[H]=get_polytopeP(lambda,F,G,Theta_HC0,p,K,A_0,A_1,A_2,B_0,B_1,B_2);
%Compute inital contraction rate rho_theta0
rho_array = [];
[A_theta0,B_theta0] = get_AB_theta(theta_bar0,A_0,A_1,A_2,B_0,B_1,B_2);
A_cl_theta0 = A_theta0 + B_theta0*K;
for i = 1:size(H,1)
    [~,rho_theta_i] = linprog(-H(i,:)*A_cl_theta0, H, ones(size(H,1),1),[],[],[],[],options);
    rho_array = [rho_array; -rho_theta_i];
end
rho_theta0 = max(rho_array);
%Compute L_B
L_array = [];
for l = 1:(2^p)
    %Compute D(x,Kx)e_tilde_l
    D_el=get_Del(Theta_HC0.V(l,:),A_1,A_2,B_1,B_2,K);
    for i = 1:size(H,1)
        [~,L_B_il] = linprog(-H(i,:)*D_el, H, ones(size(H,1),1),[],[],[],[],options);
        L_array= [L_array; -L_B_il];
    end
end
L_B = max(L_array);
%Compute c_j and c_max
c=[];
for j=1:size(F,1)
   [~,c_j]=linprog(-F(j,:)-G(j,:)*K,H,ones(size(H,1),1),[],[],[],[],options);
    c=[c;-c_j]; 
end
c_max=max(c);
% Compute bound on disturbance
d_bar_array = [];
for i = 1:(n)
    d_bar = H*W.V(i,:)';
    d_bar_array = [d_bar_array; d_bar];
end
d_bar = max(d_bar_array);
%Compute mu for RLS 1/mu>sup_{(x,u)\inZ}||D(x,u)||^2
[~,fval] = fmincon(@(x) 1/get_sqrtnormD(x,A_1,A_2,B_1,B_2),[1;1;1],[F G],ones(size(F,1),1),[],[],[],[],[],options);
mu= floor(fval);
%Check the terminal set condition
flage=check_terminalcondition(rho_theta0,eta_0,L_B,c_max,d_bar);
if flage==false
    disp("The terminal set condition is not satisfied.")
end
%Compute the Lipschitz constant on the cost function
Z=Polyhedron(L_const,l_const);
L_cost=-inf;
for i=1:size(Z.V,1)
    x=Z.V(i,1:n)';
    v=Z.V(i,n+1:end)';
    X=x+Polyhedron(H,ones(size(H,1),1));
    for j=1:size(X.V,1)
        z=X.V(j,1:n)';
        l_z=z'*Q*z+(v+K*z)'*R*(v+K*z);
        l_x=x'*Q*x+(v+K*x)'*R*(v+K*x);
        L_cost=max(L_cost,abs(l_z-l_x));
    end
end
%Steady stage cost
s_steady=1/(1-(rho_theta0+eta_0*L_B))*d_bar;
l_maxsteady=L_cost*s_steady;
end



%% Help function
function [A,B] = get_AB_theta(p,A_0,A_1,A_2,B_0,B_1,B_2)
    A = A_0 + p(1)*A_1 + p(2)*A_2;
    B = B_0 + p(1)*B_1 + p(2)*B_2;
end

function D_el=get_Del(p,A_1,A_2,B_1,B_2,K)
A = p(1)*A_1 + p(2)*A_2;
B = p(1)*B_1 + p(2)*B_2;
D_el=A+B*K;
end

function [y_sq] = get_sqrtnormD(X,A_1,A_2,B_1,B_2)
    n = length(A_1(1,:));   % state dimension
    m = length(B_1(1,:));   % input dimension
    y = [A_1*X(1:n) + B_1*X(n+1:n+m), A_2*X(1:n) + B_2*X(n+1:n+m)];
    y_sq = norm(y)^2; 
end

function flage=check_terminalcondition(rho_theta0,eta_0,L_B,c_max,d_bar)
flage=true;
if rho_theta0+eta_0*L_B+c_max*d_bar-1>0
    flage=false;
end
end