function [X_hat_OL,X_bar_OL,V_OL,S_OL,J,X,U,time,Theta_HC,t_cpu_RAMPC]=solve_RAMPC(x0,A_0,A_1,A_2,B_0,B_1,B_2,A_star,B_star,H_w,h_w,W_V,Q,R,P,K,F,G,d_bar,L_B,c_max,c,H,B_p,n,m,N,p,Theta_HC0,theta_bar0,eta_0,rho_theta0,Delta,mu,Ts,adpation)
%import Casadi
rng(1)
import casadi.*
options = optimset('Display','none',...
    'TolFun', 1e-8,...
    'MaxIter', 10000,...
    'TolConSQP', 1e-6);
%Compute the neceassary system matricis
A_cl_0=A_0+B_0*K;
A_cl_theta1=A_1+B_1*K;
A_cl_theta2=A_2+B_2*K;
%Build the constraints for y=[x_bar;x_hat;v;s] size=n*(N+1)+n*(N+1)+m*N+(N+1)
[A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_eq_theta1_hat,A_eq_theta2_hat,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho]=get_constraints(x0,A_cl_0,B_0,B_1,B_2,A_cl_theta1,A_cl_theta2,F,G,K,c,c_max,L_B,d_bar,H,B_p,n,m,N,p);
%Cost H
H_cost=getcost(Q,R,K,P,n,N);

%MPC Iteration
%initial guess
y_init=[zeros(2*(N+1)*n+N*m,1)];
for k=0:N-1
   y_init=[y_init;(1-(rho_theta0+eta_0*L_B)^k)/(1-(rho_theta0+eta_0*L_B))*d_bar];
end
y_init=[y_init;0];
%Define the requiered variables
t=0;
xmeasure = x0;
terminate=true;
Theta_HC_t=Theta_HC0;
eta_t=eta_0;
theta_bar_t=theta_bar0;
theta_hat_t=theta_bar0;
rho_theta_t=rho_theta0;
t_cpu_RAMPC=0;
time=[];
X=[];
U=[];
S=[];
J=0;
%MPC iteration loop
while t<60
    %Set membership estimation and Point estimation
    if t>0 && adpation
       %Update Theta_HC_t
       [theta_bar_t,eta_t,Theta_HC_t,Delta]=get_updatehypercube(xmeasure,x_tminus,u_tminus,theta_bar_t,eta_t,Theta_HC0,Theta_HC_t,Delta,H_w,h_w,A_0,A_1,A_2,B_0,B_1,B_2,p,B_p,options); 
       %Update theta_hat_t
       theta_hat_t=LMSpointestimate(xmeasure,x_tminus,u_tminus,Theta_HC_t,theta_hat_t,A_0,A_1,A_2,B_0,B_1,B_2,mu,p,options);
       %Update rho_theta_t
       rho_theta_t=get_updaterho(theta_bar_t,A_0,A_1,A_2,B_0,B_1,B_2,K,H,options);
    end
    
    %Optimization
    %Update the constraints with the updated parameter
    b_eq(1:n)=xmeasure;
    b_eq(n*(N+1)+1:n*(N+2))=xmeasure;
    [A_eqt,b_eqt,A_ineqt,b_ineqt]=get_currentConstraints(eta_t,rho_theta_t,theta_bar_t,theta_hat_t,A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_eq_theta1_hat,A_eq_theta2_hat,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho);
    %Optimize the quadratic programm 
        y=MX.sym('y',n*(N+1)+n*(N+1)+m*N+(N+1));%Create the optimization variable
        objective=y'*H_cost*y;%create the objective
        con=[A_eqt;A_ineqt]*y-[b_eqt;b_ineqt]; %create constraints
        lbg=[zeros(size(A_eqt,1),1);-inf*ones(size(A_ineqt,1),1)];
        ubg=[zeros(size(A_eqt,1),1);zeros(size(A_ineqt,1),1)];
        lby=[-inf*ones(n*(N+1)+n*(N+1)+m*N+(N+1),1)];
        uby=[inf*ones(n*(N+1)+n*(N+1)+m*N+(N+1),1)];
        qp = struct('x', y, 'f', objective, 'g', con);
        options_casadi = struct;
        options_casadi.print_time = 0;
        options_casadi.ipopt.print_level = 0;
        options_casadi.ipopt.hessian_constant = 'yes';
        solver = nlpsol('solver', 'ipopt', qp, options_casadi);
        tStart = cputime;
        res = solver('x0' , y_init,... % solution guess
         'lbx', lby,...           % lower bound on x
         'ubx', uby,...           % upper bound on x
         'lbg', lbg,...           % lower bound on g
         'ubg', ubg);             % upper bound on g
        y_OL=full(res.x);
        t_cpu_RAMPC = t_cpu_RAMPC+cputime - tStart;
    %Extract the open loop solution
    x_bar_OL = reshape(y_OL(1:2*(N+1))',n,[]);
    x_hat_OL=reshape(y_OL(n*(N+1)+1:2*n*(N+1)),n,[]);
    v_OL=reshape(y_OL(2*(N+1)*n+1:2*(N+1)*n+m*N),m,[]);
    s_OL=y_OL(2*(N+1)*n+m*N+1:end);
    %Store the data
        %close loop
        time=[time;t];
        X=[X,xmeasure];
        U=[U,v_OL(1:m)+K*xmeasure];
        S=[S,s_OL(1)];
        %open loop for sample set
        X_hat_OL{t+1}=x_hat_OL;
        X_bar_OL{t+1}=x_bar_OL;
        V_OL{t+1}=v_OL;
        S_OL{t+1}=s_OL;
        J_wc{t+1}=get_worstcasecosttogo(X_hat_OL{end},S_OL{end},V_OL{end},H,P,Q,R,K);
        Theta_HC{t+1}=Theta_HC_t;
    x_tminus=X(:,end);
    u_tminus=U(:,end);
    J=J+x_tminus'*Q*x_tminus+u_tminus'*R*u_tminus;
    %Simulate the uncertain system
    xmeasure=dynamic(X(:,end),U(end),A_star,B_star,W_V);
    if solver.stats.return_status~="Solve_Succeeded"
            error('Problem is infeasible')
    end
    disp(t)
    t=t+1;
end
%Compute the real time
time=(time-1)./Ts;
end


%% Help function

function [A,B]=get_systemmatrix(p,A_0,A_1,A_2,B_0,B_1,B_2)
    A = A_0 + p(1)*A_1 + p(2)*A_2;
    B = B_0 + p(1)*B_1 + p(2)*B_2;
end

function [A_eqt,b_eqt,A_ineqt,b_ineqt]=get_currentConstraints(eta,rho,theta_bar,theta_hat,A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_eq_theta1_hat,A_eq_theta2_hat,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho)
%This function builds the constraints for the current theta and rho.
    A_eqt=A_eq+A_eq_theta1_bar*theta_bar(1)+A_eq_theta2_bar*theta_bar(2)+A_eq_theta1_hat*theta_hat(1)+A_eq_theta2_hat*theta_hat(2);
    b_eqt=b_eq;
    A_ineqt=A_ineq+A_ineq_eta*eta+A_ineq_rho*rho;
    b_ineqt=b_ineq;
end

function x_tplus=dynamic(x_t,u_t,A_star,B_star,W_V)
%This function simulates the system dynamics with the disturbance
    w = W_V(1,:)';
    disturbance=w*0.5;
    x_tplus=A_star*x_t+B_star*u_t+disturbance;
end


function [A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_eq_theta1_hat,A_eq_theta2_hat,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho]=get_constraints(x0,A_cl_0,B_0,B_1,B_2,A_cl_theta1,A_cl_theta2,F,G,K,c,c_max,L_B,d_bar,H,B_p,n,m,N,p)
%This function builds the equality and the inequality constraints.
%Equality constraints
%System dynamics for x_bar
    A_eq=[eye(n),zeros(n,n*N),zeros(n,n*(N+1)),zeros(n,m*N),zeros(n,N+1)]; %inital condition
    A_eq_theta1_bar=[zeros(n),zeros(n,n*N),zeros(n,n*(N+1)),zeros(n,m*N),zeros(n,N+1)];
    A_eq_theta2_bar=[zeros(n),zeros(n,n*N),zeros(n,n*(N+1)),zeros(n,m*N),zeros(n,N+1)];
    b_eq=[x0];
    for k=0:N-1
        A_eq=[A_eq;zeros(n,n*k),A_cl_0,-eye(n),zeros(n,n*(N-k-1)),zeros(n,n*(N+1)),zeros(n,m*(k)),B_0,zeros(n,m*(N-k-1)),zeros(n,N+1)];
        A_eq_theta1_bar=[A_eq_theta1_bar;zeros(n,n*k),A_cl_theta1,zeros(n),zeros(n,n*(N-k-1)),zeros(n,n*(N+1)),zeros(n,m*k),B_1,zeros(n,m*(N-k-1)),zeros(n,N+1)];
        A_eq_theta2_bar=[A_eq_theta2_bar;zeros(n,n*k),A_cl_theta2,zeros(n),zeros(n,n*(N-k-1)),zeros(n,n*(N+1)),zeros(n,m*k),B_2,zeros(n,m*(N-k-1)),zeros(n,N+1)];    
        b_eq=[b_eq;zeros(n,1)];
    end
    %System dynamics for x_hat
    A_eq=[A_eq;zeros(n,n*(N+1)),eye(n),zeros(n,n*N),zeros(n,m*N),zeros(n,N+1)];
    b_eq=[b_eq;x0];
    A_eq_theta1_hat=[zeros(size(A_eq_theta1_bar));zeros(n),zeros(n,n*N),zeros(n,n*(N+1)),zeros(n,m*N),zeros(n,N+1)];
    A_eq_theta2_hat=[zeros(size(A_eq_theta2_bar));zeros(n),zeros(n,n*N),zeros(n,n*(N+1)),zeros(n,m*N),zeros(n,N+1)];
    for k=0:N-1
        A_eq=[A_eq;zeros(n,n*(N+1)),zeros(n,n*k),A_cl_0,-eye(n),zeros(n,n*(N-k-1)),zeros(n,m*k),B_0,zeros(n,m*(N-k-1)),zeros(n,N+1)];
        A_eq_theta1_hat=[A_eq_theta1_hat;zeros(n,n*(N+1)),zeros(n,n*k),A_cl_theta1,zeros(n),zeros(n,n*(N-k-1)),zeros(n,m*k),B_1,zeros(n,m*(N-k-1)),zeros(n,N+1)];
        A_eq_theta2_hat=[A_eq_theta2_hat;zeros(n,n*(N+1)),zeros(n,n*k),A_cl_theta2,zeros(n),zeros(n,n*(N-k-1)),zeros(n,m*k),B_2,zeros(n,m*(N-k-1)),zeros(n,N+1)];
        b_eq=[b_eq;zeros(n,1)];
    end
    %Inital condition s
    A_eq=[A_eq;zeros(1,2*n*(N+1)+N*m),1,zeros(1,N)];
    b_eq=[b_eq;0];
    A_eq_theta1_hat=[A_eq_theta1_hat;zeros(1,n*(N+1)+n*(N+1)+m*N+(N+1))];
    A_eq_theta2_hat=[A_eq_theta2_hat;zeros(1,n*(N+1)+n*(N+1)+m*N+(N+1))];
    A_eq_theta1_bar=[A_eq_theta1_bar;zeros(-size(A_eq_theta1_bar,1)+size(A_eq,1),size(A_eq_theta1_bar,2))];
    A_eq_theta2_bar=[A_eq_theta2_bar;zeros(-size(A_eq_theta2_bar,1)+size(A_eq,1),size(A_eq_theta2_bar,2))];
    %Inequality constraints
    %Terminal constraints
    A_ineq=[zeros(size(H,1),n*N),c_max*H,zeros(size(H,1),n*(N+1)),zeros(size(H,1),N*m),zeros(size(H,1),N),c_max*ones(size(H,1),1)];
    b_ineq=[ones(size(H,1),1)];
    %Constraints
    for k=0:N-1
        A_ineq=[A_ineq;zeros(size(F,1),n*k),F+G*K,zeros(size(F,1),n*(N-k)),zeros(size(F,1),n*(N+1)),zeros(size(F,1),m*k),G,zeros(size(F,1),m*(N-1-k)),zeros(size(F,1),k),c,zeros(size(F,1),N-k)];
        b_ineq=[b_ineq;ones(size(F,1),1)];
    end
    %Tube dynamic s
    A_ineq_eta=zeros(size(A_ineq));
    A_ineq_rho=zeros(size(A_ineq));
    for k=0:N-1
        for j=1:2^p  
                A_ineq_eta=[A_ineq_eta;...
              zeros(size(H,1),n*k),H*(B_p.V(j,1)*A_cl_theta1+B_p.V(j,2)*A_cl_theta2),zeros(size(H,1),n*(N-k)),zeros(size(H,1),n*(N+1)),zeros(size(H,1),m*k),H*(B_p.V(j,1)*B_1+B_p.V(j,2)*B_2),zeros(size(H,1),m*(N-1-k)),zeros(size(H,1),k),L_B*ones(size(H,1),1),zeros(size(H,1),N-k)];
        end
        b_ineq=[b_ineq;-d_bar*ones(2^p*size(H,1),1)];
        A_ineq_rho=[A_ineq_rho;...
             zeros(2^p*size(H,1),2*n*(N+1)+m*N),zeros(2^p*size(H,1),k),ones(2^p*size(H,1),1),zeros(2^p*size(H,1),N-k)];
        A_ineq=[A_ineq;...
             zeros(2^p*size(H,1),2*n*(N+1)+m*N),zeros(2^p*size(H,1),k+1),-ones(2^p*size(H,1),1),zeros(2^p*size(H,1),N-k-1)];
    end
end

function H_cost=getcost(Q,R,K,P,n,N)
%This function builds the cost function matrix H
    H_cost=[];
    %Zeros for x_bar
    H_cost=blkdiag(H_cost,zeros(n*(N+1)));
    %Cost for x_hat
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
       H_cost(2*(N+1)*n+k,(N+1)*n+k*n+1:(N+1)*n+(k+1)*n) = R*K;
       H_cost((N+1)*n+k*n+1:(N+1)*n+(k+1)*n,2*(N+1)*n+k) = K'*R;
    end
    %Zeros for s
    H_cost=blkdiag(H_cost,zeros(N+1));
    H_cost=2*H_cost; %Since we solve the problem min 0.5*y'*H*y
 end
    
function D=get_D(x,u,A_1,A_2,B_1,B_2)
    D = [A_1*x + B_1*u, A_2*x + B_2*u];
end

function [theta_bar_t,eta_t,Theta_HC_t,Delta]=get_updatehypercube(xmeasure,x_tminus,u_tminus,theta_bar_t,eta_t,Theta_HC0,Theta_HC_t,Delta,H_w,h_w,A_0,A_1,A_2,B_0,B_1,B_2,p,B_p,options)
    %This function updates the hypercube with Algorithm 1
    Theta_HC_tminus=Theta_HC_t;
    eta_tminus = eta_t;
    theta_bar_tminus = theta_bar_t;
    %Update Delta
    for i=1:length(Delta)-1
        Delta{i}=Delta{i+1};
    end
    %Compute the latest Delta
    h_delta = h_w - H_w*(xmeasure - A_0*x_tminus - B_0*u_tminus);
    H_delta = -H_w*get_D(x_tminus,u_tminus,A_1,A_2,B_1,B_2);
    Delta{end} = Polyhedron(H_delta,h_delta);
    Delta{end}.minHRep;
    %Compute Theta_M_t
    Theta_M_t = Theta_HC0; %Why
    for i = 1:length(Delta)
        Theta_M_t = Theta_M_t & Delta{i};
        Theta_M_t.minHRep;
    end
    H_theta_M_t=[Theta_M_t.A;Theta_HC_tminus.A];
    h_theta_M_t=[Theta_M_t.b;Theta_HC_tminus.b];
    Theta_M_t=Polyhedron(H_theta_M_t,h_theta_M_t);
    %Compute updated hypercube Theta_HC_t
        %Solve LPs
        for i = 1 : p
            e_i=[zeros(i-1,1);1;zeros(p-i,1)];
            [~,theta_min(i)]=linprog(e_i,H_theta_M_t,h_theta_M_t,[],[],[],[],options);
            [~,theta_max(i)]=linprog(-e_i,H_theta_M_t,h_theta_M_t,[],[],[],[],options);
            theta_max(i) = -theta_max(i);
        end
        eta_t = round(max(theta_max - theta_min),5);
        %Compute the new theta_bar_t
        theta_bar_t = 0.5*(theta_max + theta_min)';
        for i=1:p   %Projection    
           if theta_bar_t(i)< theta_bar_tminus(i)-0.5*(eta_tminus-eta_t)
                theta_bar_t(i)=theta_bar_tminus(i)-0.5*(eta_tminus-eta_t);
           elseif theta_bar_t(i)>theta_bar_tminus(i)+0.5*(eta_tminus-eta_t)
                theta_bar_t(i)=theta_bar_tminus(i)+0.5*(eta_tminus-eta_t);             
           else
             %no projection necessary  
           end          
        end
    Theta_HC_t=theta_bar_t+eta_t*B_p;
end

function theta_hat_t=LMSpointestimate(xmeasure,x_tminus,u_tminus,Theta_HC_t,theta_hat_t,A_0,A_1,A_2,B_0,B_1,B_2,mu,p,options)
%This function performs the point estimate
    theta_hat_tminus=theta_hat_t;
    [A_theta_hat_t,B_theta_hat_t]=get_systemmatrix(theta_hat_tminus,A_0,A_1,A_2,B_0,B_1,B_2);
    %predicted state
    x_hat_1t=A_theta_hat_t*x_tminus+B_theta_hat_t*u_tminus;
    %Update theta_tilde
    theta_tilde_t = theta_hat_tminus + mu*get_D(x_tminus,u_tminus,A_1,A_2,B_1,B_2)'*(xmeasure-x_hat_1t);
    [theta_hat_t,~] = fmincon(@(theta) norm(theta-theta_tilde_t),zeros(p,1),Theta_HC_t.A,Theta_HC_t.b,[],[],[],[],[],options);
end

function rho_theta_t=get_updaterho(theta_bar_t,A_0,A_1,A_2,B_0,B_1,B_2,K,H,options)
%This function updates the contraction rate for the updated theta_bar_t
    rho_array = [];
    [A_theta_bar_t,B_theta_bar_t] = get_systemmatrix(theta_bar_t,A_0,A_1,A_2,B_0,B_1,B_2);
    A_cl_theta0 = A_theta_bar_t + B_theta_bar_t*K;
    for i = 1:size(H,1)
        [~,rho_theta_i] = linprog(-H(i,:)*A_cl_theta0, H, ones(size(H,1),1),[],[],[],[],options);
        rho_array = [rho_array; -rho_theta_i];
    end
    rho_theta_t = max(rho_array);
end

function J_wc_t=get_worstcasecosttogo(X_hat_OL,S_hat_OL,V_OL,H,P,Q,R,K)
%Compute the worst-case cost-to-go
%terminal cost
X_k_t=Polyhedron(H,S_hat_OL(end)*ones(size(H,1),1)+H*X_hat_OL(:,end));
cost=0;
for i=1:size(H,1)
    cost=max([cost;X_k_t.V(i,:)*P*X_k_t.V(i,:)']);
end
cost_to_go_wc=cost;
%stage cost
for k=length(X_hat_OL)-1:-1:2
    X_k_t=Polyhedron(H,S_hat_OL(k)*ones(size(H,1),1)+H*X_hat_OL(:,k));
    cost=0;
    for i=1:size(H,1)
    cost=max([cost;X_k_t.V(i,:)*Q*X_k_t.V(i,:)'+(X_k_t.V(i,:))*K'*R*K*X_k_t.V(i,:)'+V_OL(:,k)*R*V_OL(:,k)]);
    end
    cost_to_go_wc=[cost_to_go_wc,cost_to_go_wc(end)+cost];
end
J_wc_t=flip(cost_to_go_wc);
end