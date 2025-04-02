function [H]=get_polytopeP(lambda,F,G,Theta_HC0,p,K,A_0,A_1,A_2,B_0,B_1,B_2)
options = optimset('Display','off');
%Build symmertic constrains
for i=1:size(F,2)
    F(F(:,i)~=0)=sign(F(F(:,i)~=0))*max(abs(F(F(:,i)~=0)));
end
for i=1:size(G,2)
    G(G(:,i)~=0)=sign(G(G(:,i)~=0))*max(abs(G(G(:,i)~=0)));
end
%Sift G !!!manuel %not neceassary since we control to zero
%G(5)=1/4;
%Compute the constrains for K
V=F+G*K;
L=V;
q=size(L,1);
%Compute lambda contractive set
n = 0;
fmax = -inf;
notFinished = true;    
while(notFinished)
    for j = 1:2^p
        vertices_theta=Theta_HC0.V(j,:);
        A_K = get_A_K(vertices_theta,K,A_0,A_1,A_2,B_0,B_1,B_2)/lambda;
        if n == 0
            for i = 1:q
                [~,fval] = linprog(-V(end-q+i,:)*A_K, V, ones(size(V,1),1),[],[],[],[],options);
                fmax = max(-fval-1, fmax);
            end
        elseif n > 0
            for i = 1:q*2^p
                [~,fval] = linprog(-V(end-q*2^p+i,:)*A_K, V, ones(size(V,1),1),[],[],[],[],options);
                fmax = max(-fval-1, fmax);
            end    
        end
    end
    
    if (fmax <= 0)           
        notFinished = 0;
    else
        fmax = -inf;
        Vend = [];
        for j = 1:2^p
            vertices_theta=Theta_HC0.V(j,:);
            A_K = get_A_K(vertices_theta,K,A_0,A_1,A_2,B_0,B_1,B_2)/lambda;
            Vend = [Vend; V(end-q+1:end,:)*A_K];              % To V we add A_k^i*V (another "constraint")
        end
        V = [V; Vend];
    end
    n=n+1;         % Number of iterations to find the MRPI set
end
%Compute the set and simplify
X_0 = Polyhedron(V,ones(length(V(:,1)),1));
X_0.minHRep;
X_0.minVRep;
Hx = X_0.A;
hx = X_0.b;
%Normalize constraints 
c_0 = length(Hx(:,1));
for i = 1:c_0
    H(i,:) = Hx(i,:)/hx(i);
    h(i) = 1;
end
end


%% Help function
function [A_K] = get_A_K(p,K,A_0,A_1,A_2,B_0,B_1,B_2)
    A = A_0 + p(1)*A_1 + p(2)*A_2 ;
    B = B_0 + p(1)*B_1 + p(2)*B_2 ;
    A_K = A + B*K;
end