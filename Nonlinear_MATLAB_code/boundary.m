%**************************************************************************
% Boundary Conditions
%**************************************************************************
function [K2,F1,V1,M1] = boundary(K,F,V,M,bcs)
n = size(K,1);
m = size(K,2);
q = 1;
 for j = 1:n
    if bcs(j) == 1 
        K1(:,q)= K(:,j);
        F1(q) = F(j);
        V1(q) = V(j);
        M1(q) = M(j);
        q = q + 1;
    end 
 end
 F1 = F1'; V1 = V1'; M1 = M1';
 q=1;
for j=1:n
    if bcs(j) == 1 
        K2(q,:)= K1(j,:);
        q = q + 1; 
    end 
end