function [SS_0,J_wc_0,X_bar,S,J_0]=get_Initalsolution(x0,s0,theta_bar0,B_p,eta_0,rho_theta0,L_B,d_bar,c,c_max,H,A_0,A_1,A_2,B_0,B_1,B_2,K,Q,R,P,F,G,m,n,p,L_cost,l_maxsteady,s_steady)
%compute the system matrix for theta_0
[A_cl_thetabar0,B_thetabar0] = get_AB_theta(theta_bar0,A_0,A_1,A_2,B_0,B_1,B_2);
A_cl_thetahat0=A_cl_thetabar0;
B_thetahat0=B_thetabar0;
%Simulate the dynamics
X_bar=[x0];
X_hat=[x0];
S=[s0];
V=[50.3768149142614	45.8626552601615	44.4159657269520	42.7053689031689	38.7910566037120	34.7228539889482	30.5162394242805	26.2134716705247	21.8796222120964	17.5991991995355	13.4730558543353	9.61536347171284	6.15049591728340	3.20972481169468	0.927666723099536	-0.561541362723607	-1.12833220041076	-0.652031922750973	0.509578616623437	0.516718084361636	-0.0324037289882443	-0.450780285773089	-0.640897240313527	-0.635452365059784	1.76995056610817];
k=1;
while k<200
    if k<=length(V)
        v=V(k);
    else 
        v=0;
    end
    [x_kplus_bar,x_kplus_hat,s_kplus]=dynamics(X_bar(:,k),X_hat(:,k),S(:,k),v,A_cl_thetabar0,B_thetabar0,A_cl_thetahat0,B_thetahat0,K,B_p,eta_0,rho_theta0,L_B,d_bar,H,A_1,A_2,B_1,B_2,p);
    X_bar=[X_bar,x_kplus_bar];
    X_hat=[X_hat,x_kplus_hat];
    S=[S,s_kplus];
    if check_terminalcondition(X_bar(:,k),S(:,k),s_steady)
        break;
    end
    k=k+1;
end

%Compute the worst-case cost-to-go
%terminal cost

cost_to_go_wc=0;
%stage cost
for k=length(X_bar)-1:-1:1
    x=X_bar(:,k);
    s=S(:,k);
    cost=0;
    if k<=length(V)
        v=V(k);
    else 
        v=0;
    end
    cost=x'*Q*x+(K*x+v)'*R*(K*x+v)+L_cost*s-l_maxsteady;
    cost_to_go_wc=[cost_to_go_wc,cost_to_go_wc(end)+cost];
end
%Check constraints
issatisfied=true;
v=[];
for k=1:length(X_bar)
    if k<=length(V)
        v(k)=V(k);
    else 
        v(k)=0;
    end
    
    if any((F+G*K)*X_bar(:,k)+G*v(k)+c*S(k)-1>0)
        issatisfied=false;
    end
end
if issatisfied==false
    disp("The inital solution violates the robust constraints")
end
%Construct inital sample set: 
J_wc_0=flip(cost_to_go_wc);
SS_0=[X_bar(:,1:end);S(:,1:end);v(:,1:end)];
%Compute the actual cost of the open loop trajectory
J_0=0;
for i=1:length(X_bar)
    u_k=v(i)+K*X_bar(:,i);
    J_0=J_0+X_bar(:,i)'*Q*X_bar(:,i)+u_k'*R*u_k;
end
end



%% Help function
function [A,B] = get_AB_theta(p,A_0,A_1,A_2,B_0,B_1,B_2)
    A = A_0 + p(1)*A_1 + p(2)*A_2;
    B = B_0 + p(1)*B_1 + p(2)*B_2;
end

function [x_kplus_bar,x_kplus_hat,s_kplus]=dynamics(x_k_bar,x_k_hat,s_k,v_k,A_cl_thetabar0,B_thetabar0,A_cl_thetahat0,B_thetahat0,K,B_p,eta_0,rho_theta0,L_B,d_bar,H,A_1,A_2,B_1,B_2,p)
%Input
u_k_hat=v_k+K*x_k_hat;
u_k_bar=v_k+K*x_k_bar;
%system dynamic
x_kplus_bar=A_cl_thetabar0*x_k_bar+B_thetabar0*u_k_bar;
x_kplus_hat=A_cl_thetahat0*x_k_hat+B_thetahat0*u_k_hat;
%Overapproximate the uncertainty
w_k=get_overappuncertainty(x_k_bar,u_k_bar,s_k,B_p,eta_0,rho_theta0,L_B,H,d_bar,A_1,A_2,B_1,B_2,p);
%scalar tube dynamics
s_kplus=rho_theta0*s_k+w_k;
end

function w_k=get_overappuncertainty(x_k_bar,u_k_bar,s_k,B_p,eta_0,rho_theta0,L_B,H,d_bar,A_1,A_2,B_1,B_2,p)
w_array=[];
D=get_D(x_k_bar,u_k_bar,A_1,A_2,B_1,B_2);
for l=1:2^p
    for i=1:size(H,1)
    w_array=[w_array;d_bar+eta_0*(L_B*s_k+H(i,:)*D*B_p.V(l,:)')];
    end
end
w_k=max(w_array);
end

function D=get_D(x_k_bar,u_k_bar,A_1,A_2,B_1,B_2)
D=[A_1*x_k_bar+B_1*u_k_bar,A_2*x_k_bar+B_2*u_k_bar];
end

function terminate=check_terminalcondition(x,s,s_steady)
terminate=all(abs([x;s]-[zeros(size(x,1),1);s_steady])<=10^-8);
end