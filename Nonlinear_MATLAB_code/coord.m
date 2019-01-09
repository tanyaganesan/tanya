function [cordx,cordy,cordx_node,cordy_node] = coord(dxl,dx,dxr,dz,H,Lx_l,L,Lx_r,Ly)
X1_L = [dxl/2:dxl:Lx_l-dxl/2];
X1_B = [dx/2:dx:L-dx/2];
X1_R = [dxr/2:dxr:Lx_r-dxr/2];
Y1 = [dz/2:dz:H-dz/2];

[X1L,Y1L] = meshgrid(X1_L,Y1);
[X1B,Y1B] = meshgrid(X1_B,Y1);
[X1R,Y1R] = meshgrid(X1_R,Y1);

cordx_b=(X1L(1,end)+dxl*0.5+X1B);
% Cordinate at mid point of the elements
cordx = [X1L  cordx_b (cordx_b(1,end) + dx*0.5 + X1R)];
cordy = [Y1L Y1B Y1R];

% Cordinate at node points of the elements
X1_L_n = [0:dxl:Lx_l];
X1_B_n = [0:dx:L];
X1_R_n = [0:dxr:Lx_r];
Y1_n = [0:dz:H];

[X1Ln,Y1Ln] = meshgrid(X1_L_n,Y1_n);
[X1Bn,Y1Bn] = meshgrid(X1_B_n,Y1_n);
[X1Rn,Y1Rn] = meshgrid(X1_R_n,Y1_n);

cordx_b_=(X1Ln(1,end) +  X1Bn(:,2:end-1));
cordx_node = [X1Ln  cordx_b_ (cordx_b_(1,end) + dx + X1Rn)];
cordy_node = [Y1Ln Y1Bn(:,2:end) Y1Rn(:,2:end)];
return