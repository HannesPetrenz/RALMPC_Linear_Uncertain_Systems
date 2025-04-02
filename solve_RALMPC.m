function [X,U,S,J,X_OL,V_OL,S_OL,J_wc,time,Theta_HC,t_cpu,t_solve_RALMPC]=solve_RALMPC(x0,SS,Q_func,A_0,A_1,A_2,B_0,B_1,B_2,A_star,B_star,H_w,h_w,W_V,Q,R,K,F,G,d_bar,L_B,c,H,B_p,n,m,N,p,Theta_HC0,theta_bar0,eta_0,rho_theta0,Delta,L_cost,l_maxsteady,Ts,numberitertions,adaption)
%import Casadi
import casadi.*
options = optimset('Display','none',...
    'TolFun', 1e-10,...
    'MaxIter', 1000,...
    'TolConSQP', 1e-8);
%Compute the neceassary system matricis
A_cl_0=A_0+B_0*K;
A_cl_theta1=A_1+B_1*K;
A_cl_theta2=A_2+B_2*K;
%Build the constraints for y=[x_bar;v;s;s_tilde;lamda]
%size=n*(N+1)+m*N+(N+1)+1+L
L=size(Q_func,2);
SS_init=SS;
Q_func_init=Q_func;
[A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho]=get_constraints(x0,A_cl_0,B_0,B_1,B_2,A_cl_theta1,A_cl_theta2,F,G,K,c,L_B,d_bar,H,B_p,n,m,N,p,L);
%Cost H
H_cost_base=getcost(Q,R,K,n,N,m,L_cost);
%initial guess
y_init=[zeros((N+1)*n+N*m,1)];
for k=0:N-1
   y_init=[y_init;(1-(rho_theta0+eta_0*L_B)^k)/(1-(rho_theta0+eta_0*L_B))*d_bar];
end
%Define the requiered variables
Theta_HC_t=Theta_HC0;
eta_t=eta_0;
theta_bar_t=theta_bar0;
rho_theta_t=rho_theta0;
for h=1:numberitertions
    t=0;
    time_array=[];
    X_array=[];
    U_array=[];
    S_array=[];
    J_array=0;
    xmeasure = x0;
    tStart = cputime;
    t_solve_RALMPC{h}=0;
    while t<60
        %Set membership estimation and Point estimation
        if t>0 && adaption
           %Update Theta_HC_t
           [theta_bar_t,eta_t,Theta_HC_t,Delta]=get_updatehypercube(xmeasure,x_tminus,u_tminus,theta_bar_t,eta_t,Theta_HC0,Theta_HC_t,Delta,H_w,h_w,A_0,A_1,A_2,B_0,B_1,B_2,p,B_p,options); 
           %Update rho_theta_t
           rho_theta_t=get_updaterho(theta_bar_t,A_0,A_1,A_2,B_0,B_1,B_2,K,H,options);
        end
        %Optimize 
        %Update the constraints with the updated parameter
        [A_eqt,b_eqt,A_ineqt,b_ineqt]=get_currentConstraints(SS,xmeasure,eta_t,rho_theta_t,theta_bar_t,A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho,H,L,n,N);
        %Update cost function
        [f_cost,H_cost]=get_updateCost(Q_func,n,m,N,H_cost_base,L);
        %Optimize the quadratic programm 
        if t>=0 && h>1
            y_init=Y_OL{t+1};
        elseif t>0 && h==1
            y_OL=Y_OL{end};
            y_init=[y_OL(n+1:n*(N+1),1);zeros(n,1);y_OL(n*(N+1)+m+1:n*(N+1)+m*N,1);zeros(m,1);y_OL(n*(N+1)+m*N+1:n*(N+1)+m*N+N,1);0;y_OL(n*(N+1)+m*N+N+2);y_OL(n*(N+1)+m*N+N+1+1+1:end)];
        else
           for k=0:N-1
                y_init=[y_init;(1-(rho_theta0+eta_0*L_B)^k)/(1-(rho_theta0+eta_0*L_B))*d_bar];
            end 
        end
        y_init=[y_init;zeros(n*(N+1)+m*N+(N+1)+1+L-size(y_init,1),1)];
        y_init=y_init(1:n*(N+1)+m*N+(N+1)+1+L);
        %Create optimization problem for casadi 
        y=MX.sym('y',n*(N+1)+m*N+(N+1)+1+L);%Create the optimization variable
        objective=y'*H_cost*y+f_cost'*y;%create the objective
        con=[A_eqt;A_ineqt]*y-[b_eqt;b_ineqt]; %create constraints
        lbg=[zeros(size(A_eqt,1),1);-inf*ones(size(A_ineqt,1),1)];
        ubg=[zeros(size(A_eqt,1),1);zeros(size(A_ineqt,1),1)];
        lby=[-inf*ones(n*(N+1)+m*N+(N+1)+1,1);zeros(L,1)];
        uby=[inf*ones(n*(N+1)+m*N+(N+1)+1,1);ones(L,1)];
        qp = struct('x', y, 'f', objective, 'g', con);
        options_casadi = struct;
        options_casadi.print_time = 0;
        options_casadi.ipopt.print_level = 0;
        %options_casadi.ipopt.hessian_constant = 'yes';
        tStart_solve = cputime;
        solver = nlpsol('solver', 'ipopt', qp, options_casadi);
        res = solver('x0' , y_init,... % solution guess
         'lbx', lby,...           % lower bound on x
         'ubx', uby,...           % upper bound on x
         'lbg', lbg,...           % lower bound on g
         'ubg', ubg);             % upper bound on g
        t_solve_RALMPC{h}=t_solve_RALMPC{h}+cputime - tStart_solve;
        y_OL=full(res.x);
        %Extract the open loop solution
        x_OL = reshape(y_OL(1:n*(N+1))',n,[]);
        v_OL=reshape(y_OL((N+1)*n+1:(N+1)*n+m*N),m,[]);
        s_OL=y_OL((N+1)*n+m*N+1:(N+1)*n+m*N+N+1);
        lambda=y_OL((N+1)*n+m*N+N+3:end);
        Q_wc=f_cost'*y_OL; 
        %Store the data
        Y_OL{t+1}=y_OL;
        time_array=[time_array;t];
        X_array=[X_array,xmeasure];
        U_array=[U_array,v_OL(1,:)+K*xmeasure];
        S_array=[S_array,s_OL(1,:)];
        %open loop for sample set
        X_OL{h}{t+1}=x_OL;
        V_OL{h}{t+1}=v_OL;
        S_OL{h}{t+1}=s_OL;
        J_wc{h}{t+1}=get_worstcasecosttogo(x_OL,s_OL,v_OL,Q_wc,Q,R,K,L_cost,l_maxsteady);
        Theta_HC{h}{t+1}=Theta_HC_t;
        x_tminus=xmeasure;
        u_tminus=v_OL(:,1)+K*xmeasure;
        %Simulate the uncertain system
        xmeasure=dynamic(x_tminus,u_tminus,A_star,B_star,W_V);
        %Compute cost function
        J_array=J_array+x_tminus'*Q*x_tminus+u_tminus'*R*u_tminus;
        %Update time
        t=t+1;
        if solver.stats.return_status~="Solve_Succeeded"
            error('Problem is infeasible')
        end
    end
    %Compute the real time
    time_array=(time_array-1)./Ts;
    %close loop
    time{h}=time_array;
    X{h}=X_array;
    U{h}=U_array;
    S{h}=S_array;
    J{h}=J_array;
    %Update the sample set and the function Q
    [SS,Q_func,L]=get_updateSampleSet(SS,Q_func,X_OL{h},J_wc{h},V_OL{h},S_OL{h},Q_func_init,SS_init);
    t_cpu{h} = cputime - tStart;
end
end

%% Help functions
function D=get_D(x,u,A_1,A_2,B_1,B_2)
    D = [A_1*x + B_1*u, A_2*x + B_2*u];
end

function [A,B]=get_systemmatrix(p,A_0,A_1,A_2,B_0,B_1,B_2)
    A = A_0 + p(1)*A_1 + p(2)*A_2;
    B = B_0 + p(1)*B_1 + p(2)*B_2;
end


function x_tplus=dynamic(x_t,u_t,A_star,B_star,W_V)
%This function simulates the system dynamics with the disturbance
    w = W_V(1,:)';
    disturbance=w*0.5;
    x_tplus=A_star*x_t+B_star*u_t+disturbance;
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

function J_wc_t=get_worstcasecosttogo(X_OL,S_OL,V_OL,Q_wc,Q,R,K,L_cost,l_maxsteady)
%Compute the worst-case cost-to-go
%terminal cost
    cost_to_go_wc=Q_wc;%Set the terminal worst case cost to go. 
    %stage cost
    for k=length(X_OL)-1:-1:1
        x_k_t=X_OL(:,k);
        s_k_t=S_OL(k);
        v_k_t=V_OL(:,k);
        cost=0;
        cost=x_k_t'*Q*x_k_t+(K*x_k_t+v_k_t)'*R*(K*x_k_t+v_k_t)+L_cost*s_k_t-l_maxsteady;
        cost_to_go_wc=[cost_to_go_wc,cost_to_go_wc(end)+cost];
    end
    J_wc_t=flip(cost_to_go_wc);
end


function [A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho]=get_constraints(x0,A_cl_0,B_0,B_1,B_2,A_cl_theta1,A_cl_theta2,F,G,K,c,L_B,d_bar,H,B_p,n,m,N,p,L)
%This function builds the equality and the inequality constraints.
%Equality constraints
%System dynamics for x_bar
    A_eq=[eye(n),zeros(n,n*N),zeros(n,m*N),zeros(n,N+1),zeros(n,1)]; %inital condition
    A_eq_theta1_bar=[zeros(n),zeros(n,n*N),zeros(n,m*N),zeros(n,N+1),zeros(n,1)];
    A_eq_theta2_bar=[zeros(n),zeros(n,n*N),zeros(n,m*N),zeros(n,N+1),zeros(n,1)];
    b_eq=[x0];
    for k=0:N-1
        A_eq=[A_eq;zeros(n,n*k),A_cl_0,-eye(n),zeros(n,n*(N-k-1)),zeros(n,m*(k)),B_0,zeros(n,m*(N-k-1)),zeros(n,N+1),zeros(n,1)];
        A_eq_theta1_bar=[A_eq_theta1_bar;zeros(n,n*k),A_cl_theta1,zeros(n),zeros(n,n*(N-k-1)),zeros(n,m*k),B_1,zeros(n,m*(N-k-1)),zeros(n,N+1),zeros(n,1)];
        A_eq_theta2_bar=[A_eq_theta2_bar;zeros(n,n*k),A_cl_theta2,zeros(n),zeros(n,n*(N-k-1)),zeros(n,m*k),B_2,zeros(n,m*(N-k-1)),zeros(n,N+1),zeros(n,1)];    
        b_eq=[b_eq;zeros(n,1)];
    end
    %Inital condition s
    A_eq=[A_eq;zeros(1,n*(N+1)+N*m),1,zeros(1,N),zeros(1,1)];
    b_eq=[b_eq;0];
    A_eq_theta1_bar=[A_eq_theta1_bar;zeros(-size(A_eq_theta1_bar,1)+size(A_eq,1),size(A_eq_theta1_bar,2))];
    A_eq_theta2_bar=[A_eq_theta2_bar;zeros(-size(A_eq_theta2_bar,1)+size(A_eq,1),size(A_eq_theta2_bar,2))];
    %Inequality constraints
    %Terminal constraints
    A_ineq=[zeros(size(H,1),n*N),H,zeros(size(H,1),N*m),zeros(size(H,1),N+1),-ones(size(H,1),1)];
    b_ineq=zeros(size(H,1),1);
    A_ineq=[A_ineq;zeros(1,n*(N+1)),zeros(1,N*m),zeros(1,N),ones(1,1),ones(1,1)];
    b_ineq=[b_ineq;0];
    b_ineq=[b_ineq;1];%lambda<=1
    b_ineq=[b_ineq;zeros(L,1)];%lambda_i>0
    A_ineq=[A_ineq;zeros(L+1,n*(N+1)),zeros(L+1,N*m),zeros(L+1,N+2)];
    %Constraints
    for k=0:N-1
        A_ineq=[A_ineq;zeros(size(F,1),n*k),F+G*K,zeros(size(F,1),n*(N-k)),zeros(size(F,1),m*k),G,zeros(size(F,1),m*(N-1-k)),zeros(size(F,1),k),c,zeros(size(F,1),N-k),zeros(size(F,1),1)];
        b_ineq=[b_ineq;ones(size(F,1),1)];
    end
    %Tube dynamic s
    A_ineq_eta=zeros(size(A_ineq));
    A_ineq_rho=zeros(size(A_ineq));
    for k=0:N-1
        for j=1:2^p  
                A_ineq_eta=[A_ineq_eta;...
              zeros(size(H,1),n*k),H*(B_p.V(j,1)*A_cl_theta1+B_p.V(j,2)*A_cl_theta2),zeros(size(H,1),n*(N-k)),zeros(size(H,1),m*k),H*(B_p.V(j,1)*B_1+B_p.V(j,2)*B_2),zeros(size(H,1),m*(N-1-k)),zeros(size(H,1),k),L_B*ones(size(H,1),1),zeros(size(H,1),N-k),zeros(size(H,1),1)];
        end
        b_ineq=[b_ineq;-d_bar*ones(2^p*size(H,1),1)];
        A_ineq_rho=[A_ineq_rho;...
             zeros(2^p*size(H,1),n*(N+1)+m*N),zeros(2^p*size(H,1),k),ones(2^p*size(H,1),1),zeros(2^p*size(H,1),N-k),zeros(2^p*size(H,1),1)];
        A_ineq=[A_ineq;...
             zeros(2^p*size(H,1),n*(N+1)+m*N),zeros(2^p*size(H,1),k+1),-ones(2^p*size(H,1),1),zeros(2^p*size(H,1),N-k-1),zeros(2^p*size(H,1),1)];
    end   
end

function H_cost=getcost(Q,R,K,n,N,m,L_cost)
%This function builds the cost function matrix H
    H_cost=[];
    %Cost for x
    for k=0:N-1
       H_cost=blkdiag(H_cost,Q+K'*R*K); 
    end
    %Terminal cost
    H_cost=blkdiag(H_cost,zeros(size(Q)));
    %Cost for v
    for k=0:N-1
       H_cost=blkdiag(H_cost,R); 
    end
    %Cost for s
    for k=0:N-1
       H_cost=blkdiag(H_cost,L_cost); 
    end
    %s_N
    H_cost=blkdiag(H_cost,zeros(1));
    %s_tilde
    H_cost=blkdiag(H_cost,zeros(1));
    %Add cross-terms
    for k=0:N-1
        H_cost((N+1)*n+k+m,k*n+1:(k+1)*n) = R*K;
        H_cost(k*n+1:(k+1)*n,(N+1)*n+k+1) = K'*R;
    end
    H_cost=2*H_cost; %Since we solve the problem min 0.5*y'*H*y
end
function [f_cost,H_cost]=get_updateCost(Q_func,n,m,N,H_cost_base,L)
f_cost=[zeros(n*(N+1)+m*N+(N+1)+1,1);Q_func'];
H_cost=blkdiag(H_cost_base,zeros(L));
end

function [SS,Q_func,L]=get_updateSampleSet(SS,Q_func,X_bar,J_wc,V_OL,S_OL,Q_func_init,SS_init)
    %This function updates the Sample set and Q_fun
    for i=1:length(X_bar)
        SS=[SS,[X_bar{i}(:,1:end-1);S_OL{i}(1:end-1)';V_OL{i}(:,1:end)]];
        Q_func=[Q_func,J_wc{i}(1:end)];
    end
    %Convex Hull
    idx=convhull(SS(1,:),SS(2,:));
    SS=SS(:,idx);
    Q_func=Q_func(idx);
    %Add initial solution for feasibility
    [C,ia,ic]=unique([SS,SS_init]','stable',"rows");
    SS=C';
    Q_func=[Q_func,Q_func_init]';
    Q_func=Q_func(ia)';
    %Required length lambda
    L=size(Q_func,2);
end

function [A_eqt,b_eqt,A_ineqt,b_ineqt]=get_currentConstraints(SS,x0,eta,rho,theta_bar,A_eq,b_eq,A_eq_theta1_bar,A_eq_theta2_bar,A_ineq,b_ineq,A_ineq_eta,A_ineq_rho,H,L,n,N)
%This function builds the constraints for the current theta and rho.
    b_eq(1:n)=x0;
    A_eqt=A_eq+A_eq_theta1_bar*theta_bar(1)+A_eq_theta2_bar*theta_bar(2);
    b_eqt=b_eq;
    A_ineqt=A_ineq+A_ineq_eta*eta+A_ineq_rho*rho;
    b_ineqt=b_ineq;
    %Update Lambda part
    A_ineq_lamda=-H*SS(1:2,:);
    A_ineq_lamda=[A_ineq_lamda;-SS(3,:)];
    A_ineq_lamda=[A_ineq_lamda;ones(1,L);-eye(L)];
    A_ineq_lamda=[A_ineq_lamda;zeros(size(A_ineq,1)-size(A_ineq_lamda,1),L)];

    A_eq_lambda=zeros(size(b_eq,1),L);

    A_ineqt=[A_ineqt,A_ineq_lamda];
    A_eqt=[A_eqt,A_eq_lambda];
end

