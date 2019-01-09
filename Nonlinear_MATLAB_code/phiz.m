function [FI_] = phiz(d_soil_L,disV,d_soil_R,b,H,SpE,nu_,leL,leb,leR,nn1,nn2,nn3,CorX,CorY,CorX_n,CorY_n)
dd_soil_L = diff(d_soil_L)/leL;
dd_beam = diff(disV)/leb;
dd_soil_R = diff(d_soil_R)/leR;
%********************
A1 = b*(1-nu_)/((1+nu_)*(1-2*nu_));
B1 = b/(2*(1+nu_));
% Es values at grid point
SpE_n = zeros(size(SpE)+1);
SpE_n(1:end-1,1:end-1) = SpE;
SpE_n(1:end-1,end) = SpE(:,end);
SpE_n(end,:) = SpE_n(end-1,:);
% *************************************************
for J = 1:size(SpE_n,1);
m(J,1) = A1*( 0.5*(SpE_n(J,1)*d_soil_L(1)^2*leL ...
+ SpE_n(J,end)*d_soil_R(end)^2*leR) + ...   
sum(SpE_n(J,2:nn1-1)'.*d_soil_L(2:end-1).^2*leL) + ...
sum(SpE_n(J,nn1:nn1+nn2-1)'.*disV.^2*leb) + ...
sum(SpE_n(J,nn1+nn2:end-1)'.*d_soil_R(2:end-1).^2*leR) );
% 
n(J,1) = B1*( 0.5*(SpE_n(J,1)*dd_soil_L(1)^2*leL + SpE_n(J,end-1)*dd_soil_R(end)^2*leR) + ...   
sum(SpE_n(J,2:nn1-1)'.*dd_soil_L(2:end).^2*leL) + sum(SpE_n(J,nn1:nn1+nn2-2)'.*dd_beam.^2*leb) + ...
sum(SpE_n(J,nn1+nn2:end-1)'.*dd_soil_R(1:end-1).^2*leR) );
end
% 
% 
N = size(SpE_n,1);
h = H/(N-1);
for j = [1:1];
    FI_z(1,j)   = 2*m(j+1)/h^2 + n(j+1);
    FI_z(1,j+1) = -m(j+1)/h^2 - 1/(4*h^2)*(m(j+2)-m(j));
    FI_z(1,j+2) = 0;
end
% 
for j = [2:N-3];
    FI_z(j,j-1) = -m(j+1)/h^2 + 1/(4*h^2)*(m(j+2)-m(j));
    FI_z(j,j)   = 2*m(j+1)/h^2 + n(j+1);
    FI_z(j,j+1) = -m(j+1)/h^2 - 1/(4*h^2)*(m(j+2)-m(j));
end 
% 
for j = [N-2:N-2]
    FI_z(j,j-1) = -m(j+1)/h^2 + 1/(4*h^2)*(m(j+2)-m(j+1));
    FI_z(j,j) = 2*m(j+1)/h^2 + n(j+1);
end
F = [m(2)/h^2 - 1/(4*h^2)*(m(3)-m(1));zeros((N-3),1)];
FI_ = inv(FI_z)*F;
FI_ = [1;FI_;0];
return