function [gm_oct,tao_xz,sig_mn,eZZ,eXZ,ev] = strss_strain(SpE,CORDY,Fz,dFzdz,n,m,nesL,neB,nesR, ...
                                             d_,disV,dxl,dx,dxr,nu_,gamma_z)
% 
% Calculation of representative shear strain
n1 = size(d_,1); m1 = size(dFzdz,1);
exx = zeros(m1,n1);
ezz = -dFzdz* d_';
% 
dW_dx = [diff(d_(1:nesL+1))/dxl; diff(disV)/dx; diff(d_(nesL+neB+1:end))/dxr];
exz = -1/2*Fz*dW_dx';
% 
gm_oct = 2/3*(2*ezz(1:n,1:m).^2 + 6*exz(1:n,1:m).^2).^0.5;
ev = ezz(1:n,1:m); % Volumetric strain (1D)
eZZ = ezz(1:n,1:m);
eXZ = exz(1:n,1:m);
% 
geo_str = CORDY*gamma_z + gamma_z*2; % Geostatic Stress

sig_xx = SpE/((1+nu_)*(1-2*nu_))*nu_/(1-nu_).*eZZ;
sig_zz = SpE*nu_/(1-nu_).*eZZ + geo_str;
sig_yy = nu_*(sig_xx + sig_zz);
tao_xz = SpE*(1-2*nu_)/(2*(1-nu_)).*eXZ;
sig_mn = 1/3*(abs(sig_xx) + abs(sig_yy) + abs(sig_zz));
return