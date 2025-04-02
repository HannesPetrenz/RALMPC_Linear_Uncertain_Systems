function [P,K]=get_P_K(rho,Q,R,A_0,A_1,A_2,B_0,B_1,B_2,Theta_HC0,m,n,p)

%Define the optimization variables
Y = sdpvar(m,n); 
X = sdpvar(n);
%Vertices theta_j
p_minmax=zeros(p,2);
for i=1:p
p_minmax(i,:)=[min(Theta_HC0.V(:,i)),max(Theta_HC0.V(:,i))];
end
theta_j=table2array(combinations(p_minmax(1,:),p_minmax(2,:))); %extend to more parameter
%Build the constrains
con=[];
for i=1:length(theta_j)
    [A_theta,B_theta]=get_ABtheta(theta_j(i,:),A_0,A_1,A_2,B_0,B_1,B_2); %Compute A_thetaj and B_thetaj
    ineq = [X, (A_theta*X+B_theta*Y)', sqrtm(Q)*X, (sqrtm(R)*Y)';...     %Equation (33b)
                A_theta*X+B_theta*Y, X, zeros(n,n+m);... 
                X*sqrtm(Q), zeros(n), eye(n), zeros(n,m);...
                sqrtm(R)*Y, zeros(m,2*n), eye(m)];
    con = [con;ineq>=0];
    ineq = [rho*X, (A_theta*X+B_theta*Y)';...                             %Equation (33c)
                A_theta*X+B_theta*Y, rho*X];
    con = [con;ineq>=0];
end
%Solve the optimization problem
optimize(con,-log(det(X))); %min_{X,Y} -logdet(X) subject to con Yalmip
%Compute P and K
Y = value(Y);
X = value(X);

P=inv(X);
K=Y*P;
%Check the solution
flage=check(P,K,Q,R,theta_j,A_0,A_1,A_2,B_0,B_1,B_2);
if flage==false
disp("Solution P and K do not satisfy Assumption 2")
end
end
%% Help functions
function [A_theta,B_theta]=get_ABtheta(theta_j,A_0,A_1,A_2,B_0,B_1,B_2)
A_theta = A_0 + theta_j(1)*A_1 + theta_j(2)*A_2 ;
B_theta = B_0 + theta_j(1)*B_1 + theta_j(2)*B_2 ;
end

function flage=check(P,K,Q,R,theta_j,A_0,A_1,A_2,B_0,B_1,B_2)
flage=true;
for i=1:length(theta_j)
    [A_theta,B_theta]=get_ABtheta(theta_j(i,:),A_0,A_1,A_2,B_0,B_1,B_2);
    A_cl_theta=A_theta+B_theta*K;
    cond=A_cl_theta'*P*A_cl_theta+Q+K'*R*K-P;
    if any(round(eig(cond),6)>0)
            flage=false;
    end
end
end