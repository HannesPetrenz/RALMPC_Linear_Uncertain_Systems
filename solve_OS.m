function [X_bar_OL,V_OL,S_OL,J,X,U,time]=solve_OS(x0,A_star,B_star,Q,R,P,K,F,G,d_bar,c_max,c,H,n,m,N,W_V,Ts)
%import Casadi
rng(1)
import casadi.*
options = optimset('Display','none',...
    'TolFun', 1e-8,...
    'MaxIter', 10000,...
    'TolConSQP', 1e-6);
%Compute the neceassary system matricis
A_cl_0=A_star+B_star*K;
rho=get_updaterho(A_cl_0,B_star,H,options);
%Build the constraints for y=[x_bar;v;s] size=n*(N+1)+n*(N+1)+m*N+(N+1)
[A_eq,b_eq,A_ineq,b_ineq]=get_constraints(x0,A_cl_0,B_star,F,G,K,c,c_max,d_bar,H,n,m,N,rho);
%Cost H
H_cost=getcost(Q,R,K,P,n,m,N);

%MPC Iteration
%initial guess
y_init=[zeros((N+1)*n+N*m,1)];
for k=0:N-1
   y_init=[y_init;d_bar];
end
y_init=[y_init;0];
%Define the requiered variables
t=0;
xmeasure = x0;

time=[];
X=[];
U=[];
S=[];
J=0;
%MPC iteration loop
while t<60
    %Optimization
    %Update the constraints with the updated parameter
    b_eq(1:n)=xmeasure;
        y=MX.sym('y',n*(N+1)+m*N+(N+1));%Create the optimization variable
        objective=y'*H_cost*y;%create the objective
        con=[A_eq;A_ineq]*y-[b_eq;b_ineq]; %create constraints
        lbg=[zeros(size(A_eq,1),1);-inf*ones(size(A_ineq,1),1)];
        ubg=[zeros(size(A_eq,1),1);zeros(size(A_ineq,1),1)];
        lby=[-inf*ones(n*(N+1)+m*N+(N+1),1)];
        uby=[inf*ones(n*(N+1)+m*N+(N+1),1)];
        qp = struct('x', y, 'f', objective, 'g', con);
        options_casadi = struct;
        options_casadi.print_time = 0;
        options_casadi.ipopt.print_level = 0;
        %options_casadi.ipopt.hessian_constant = 'yes';
        solver = nlpsol('solver', 'ipopt', qp, options_casadi);
        res = solver('x0' , y_init,... % solution guess
         'lbx', lby,...           % lower bound on x
         'ubx', uby,...           % upper bound on x
         'lbg', lbg,...           % lower bound on g
         'ubg', ubg);             % upper bound on g
        y_OL=full(res.x);
    %Extract the open loop solution
    x_bar_OL = reshape(y_OL(1:2*(N+1))',n,[]);
    v_OL=reshape(y_OL((N+1)*n+1:(N+1)*n+m*N),m,[]);
    s_OL=y_OL((N+1)*n+m*N+1:end);
    %Store the data
        %close loop
        time=[time;t];
        X=[X,xmeasure];
        U=[U,v_OL(1:m)+K*xmeasure];
        S=[S,s_OL(1)];
        %open loop for sample set
        X_bar_OL{t+1}=x_bar_OL;
        V_OL{t+1}=v_OL;
        S_OL{t+1}=s_OL;
    x_tminus=X(:,end);
    u_tminus=U(:,end);
    J=J+x_tminus'*Q*x_tminus+u_tminus'*R*u_tminus;
    %Simulate the uncertain system
    xmeasure=dynamic(X(:,end),U(:,end),A_star,B_star,W_V);
    if solver.stats.return_status~="Solve_Succeeded"
            error('Problem is infeasible')
    end
    t=t+1;
end
%Compute the real time
time=(time-1)./Ts;
end


%% Help function

function x_tplus=dynamic(x_t,u_t,A_star,B_star,W_V)
%This function simulates the system dynamics with the disturbance
    w = W_V(1,:)';
        disturbance=w*0.5;
    x_tplus=A_star*x_t+B_star*u_t+disturbance;
end


function [A_eq,b_eq,A_ineq,b_ineq]=get_constraints(x0,A_cl_0,B_0,F,G,K,c,c_max,d_bar,H,n,m,N,rho)
%This function builds the equality and the inequality constraints.
%Equality constraints
%System dynamics for x_bar
    A_eq=[eye(n),zeros(n,n*N),zeros(n,m*N),zeros(n,N+1)]; %inital condition
    b_eq=[x0];
    for k=0:N-1
        A_eq=[A_eq;zeros(n,n*k),A_cl_0,-eye(n),zeros(n,n*(N-k-1)),zeros(n,m*(k)),B_0,zeros(n,m*(N-k-1)),zeros(n,N+1)];
        b_eq=[b_eq;zeros(n,1)];
    end
    
    %Inital condition s
    A_eq=[A_eq;zeros(1,n*(N+1)+N*m),1,zeros(1,N)];
    b_eq=[b_eq;0];    %Inequality constraints
    %Terminal constraints
    A_ineq=[zeros(size(H,1),n*N),c_max*H,zeros(size(H,1),N*m),zeros(size(H,1),N),c_max*ones(size(H,1),1)];
    b_ineq=[ones(size(H,1),1)];
    %Constraints
    for k=0:N-1
        A_ineq=[A_ineq;zeros(size(F,1),n*k),F+G*K,zeros(size(F,1),n*(N-k)),zeros(size(F,1),m*k),G,zeros(size(F,1),m*(N-1-k)),zeros(size(F,1),k),c,zeros(size(F,1),N-k)];
        b_ineq=[b_ineq;ones(size(F,1),1)];
    end
    %Tube dynamic s
    for k=0:N-1
        b_ineq=[b_ineq;-d_bar];
        A_ineq=[A_ineq;...
             zeros(1,n*(N+1)+m*N),zeros(1,k),rho,-1,zeros(1,N-k-1)];
    end
end

function H_cost=getcost(Q,R,K,P,n,m,N)
%This function builds the cost function matrix H
    H_cost=[];
    %Cost for x_bar
    for k=0:N-1
       H_cost=blkdiag(H_cost,Q+K'*R*K); 
    end
    %Terminal cost
    H_cost=blkdiag(H_cost,P);
    %Cost for v
    for k=0:N-1
       H_cost=blkdiag(H_cost,R); 
    end
    %Add cross-terms
    for k=0:N-1
        H_cost((N+1)*n+k+m,k*n+1:(k+1)*n) = R*K;
        H_cost(k*n+1:(k+1)*n,(N+1)*n+k+1) = K'*R;
    end
    %Zeros for s
    H_cost=blkdiag(H_cost,zeros(N+1));
    H_cost=2*H_cost; %Since we solve the problem min 0.5*y'*H*y
 end
    
function rho_theta_t=get_updaterho(A_cl,B_0,H,options)
%This function updates the contraction rate for the updated theta_bar_t
    rho_array = [];
    for i = 1:size(H,1)
        [~,rho_theta_i] = linprog(-H(i,:)*A_cl, H, ones(size(H,1),1),[],[],[],[],options);
        rho_array = [rho_array; -rho_theta_i];
    end
    rho_theta_t = max(rho_array);
end

