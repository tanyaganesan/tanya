function [d_,disV,BM,SF,nn2,k_,t_,z,Fi_z_,dFi_dz] = ...
    analysis(nesL,neB,nesR,m,n,SpE,L,Lx,Ly, ...
    nu_,H,b,Eb,Ib,Q,V,Mo,l_pos,dxl,dx,dxr,dz,Lx_l,Lx_r)
[CorX,CorY,CorXn,CorYn] = coord(dxl,dx,dxr,dz,H,Lx_l,L,Lx_r,Ly);
% [SpE] = RandomField(nesL,neB,nesR,m,n,Es,0.5,0.5,L,Ly,CorX,CorY);
% [Sp_nu] = RandomField(nesL,neB,nesR,m,n,Mean_nu,COV_nu,Corr_Dist_nu,Lx,Ly,CorX,CorY);
% SpE = ones(n,m)*Es;
% Sp_nu = ones(n,m)*nu_;
z_mid = CorY(:,1); 
dz = H/n;
z = [0:dz:H]';
Fi_z = sinh(4*(1-z/H))/sinh(4);
dFi_z = diff(Fi_z);
G_ = b/(4*(1 + nu_));
Lamda_ = b*(1-nu_)/((1+nu_)*(1-2*nu_));
%**************************************************************************
nn1 = nesL + 1;    % Number of Nodes : soil present at left side of beam
nn2 = neB  + 1;    % Number of Nodes : beam
nn3 = nesR + 1;    % Number of Nodes : soil present at right side of beam
%**************************************************************************
lesL = (Lx-L)/2; lesR = lesL;
leL = lesL/nesL;   % Length of soil element (left)
leb = L/neB;       % Length of beam element 
leR = lesR/nesR;   % Length of soil element (right)
%**************************************************************************
element1 = [1:1:(nn1-1); 2:1:nn1]';                                 % Element Connectivity: soil (left)
element2 = [1:2:(2*nn2-3); 2:2:2*nn2-2 ; 3:2:2*nn2-1 ; 4:2:2*nn2]'; % Element Connectivity: beam
element3 = [1:1:(nn3-1); 2:1:nn3]';                                 % Element Connectivity: soil (right)
%**************************************************************************
F1 = zeros(nn2*2,1);                          % Initial global force vector for distributed load
% ************************************************
% Position of Concentrated Loads
% ************************************************
n_pos(1) = 2*round(l_pos(1)/leb) + 1;
n_pos(2) = 2*round(l_pos(2)/leb) + 1;
n_pos(3) = 2*round(l_pos(3)/leb) + 1;
% 
for II = [1:1:2*nn2];
    if II == n_pos(1) | II == n_pos(2) | II == n_pos(3)
        F2_(II,1) = 1;
    else
        F2_(II,1) = 0;
    end
end
F2 = F2_*V;                  % force vector for Concentrated Load on beam]
% ************************************************
M1 = [zeros(neB,1);0;Mo;zeros(2*nn2-neB-2,1)];   % force vector Concentrated Moment on beam 
for e4 = [1:neB];
    [qb,qs]=load_vector(leb,0.0,Q);
    ss4 = element2(e4,:);
    F1(ss4,1)=F1(ss4,1)+qb;                      % Global force vector for distributed load
end
F = [zeros(nn1-1,1);F1;zeros(nn3-1,1)];
S = [zeros(nn1-1,1);F2;zeros(nn3-1,1)];
M = [zeros(nn1-1,1);M1;zeros(nn3-1,1)];
%**************************************************************************
% Boundary conditions (0):fixed, (1):free
%**************************************************************************
%     |----soil(L)-----|---------beam----------|-----soil(R)--------|
bcs =  [1;ones(nn1-2,1);1;1;ones(2*nn2-4,1);1;1;ones(nn3-2,1);1];  %|
%     |----------------|-----------------------|--------------------|
%**************************************************************************                  
%**************************************************************************
for ITR = 1:20
Fi_z_(:,ITR) = Fi_z;
dFi_dz(:,ITR) = dFi_z/dz;
%**************************************************************************
for J = 1:m
    Fi_z_m = interp1(z,Fi_z,z_mid);
    t(J) = G_* dz * ((0.5*(SpE(1,J)*Fi_z_m(1)^2 + SpE(end,J)*Fi_z_m(end)^2) + ...
        sum(SpE(2:end-1,J).*Fi_z_m(2:end-1).^2)));
    k(J) = Lamda_*dz*( SpE(1,J)*(dFi_z(1)/dz)^2 +  SpE(end,J)*(dFi_z(end)/dz)^2  + ...
        sum(SpE(2:end-1,J).*(dFi_z(2:end-1)/dz).^2));
end
t_(:,ITR) = t(nesL:nesL+neB)';
k_(:,ITR) = k(nesL:nesL+neB)';
%**************************************************************************
K1 = zeros(nn1,nn1);                        % Initial Global Stiffness Matrix : soil (left)
K2 = zeros(nn2*2,nn2*2);                    % Initial Global Stiffness Matrix : beam 
K3 = zeros(nn3,nn3);                        % Initial Global Stiffness Matrix : soil (right)
K  = zeros(nn1+2*nn2+nn3-2,nn1+2*nn2+nn3-2);% Initial Global Stiffness Matrix : soil(left) + beam + soil (right)
%**************************************************************************
% Assembling of Stiffness Matrix
%**************************************************************************
for e1 = [1:nesL];
    [kesL1,kesL2,kesR1,kesR2] = soil_sitff(leL,t(e1),k(e1),0,0);
    ss1=element1(e1,:);
    K1(ss1,ss1)=K1(ss1,ss1)+(kesL1+kesL2);
end
for e2=[1:neB];
    [kbe1,kbe2,kbe3] = beam_sitff(leb,Eb,Ib,t(nesL+e2),k(nesL+e2));
    ss2=element2(e2,:);
    K2(ss2,ss2)=K2(ss2,ss2)+(kbe1+kbe2+kbe3);
end
for e3=[1:nesR];
    [kesL1,kesL2,kesR1,kesR2] = soil_sitff(leR,0,0,t(nesL+neB+e3),k(nesL+neB+e3));
    ss3=element3(e3,:);
    K3(ss3,ss3)=K3(ss3,ss3)+(kesR1+kesR2);
end
%**************************************************************************
% Assembling of soil(left) + beam + soil (right) stiffness matrix
%**************************************************************************
K(1:nn1,1:nn1) = K(1:nn1,1:nn1)+K1;
K(nn1:2*nn2+nn1-1,nn1:2*nn2+nn1-1) = K(nn1:2*nn2+nn1-1,nn1:2*nn2+nn1-1)+K2;
K(2*nn2+nn1-2:2:2*nn2+nn1,2*nn2+nn1-2:2:2*nn2+nn1) = ...
    K(2*nn2+nn1-2:2:2*nn2+nn1,2*nn2+nn1-2:2:2*nn2+nn1)+K3(1:2,1:2);
K3(2,2)=0;
K(2*nn2+nn1:nn1+2*nn2+nn3-2,2*nn2+nn1:nn1+2*nn2+nn3-2) = ...
    K(2*nn2+nn1:nn1+2*nn2+nn3-2,2*nn2+nn1:nn1+2*nn2+nn3-2)+K3(2:end,2:end);
%**************************************************************************
[K_,F_,S_,M_] = boundary(K,F,S,M,bcs); % Stiffness matrix & force vector after applying EBCs
u = inv(K_)*(F_+ S_+ M_);              % Global displacement vector for soil &  beam
JJ = 1;
for I = [1:nn1+2*nn2+nn3-2];
    if bcs(I) == 0
        dis(I) = 0;
    else
        dis(I) = u(JJ);
        JJ = JJ+1;
    end
end
%**************************************************************************
% Calculation of beam bending moment & shear force
%**************************************************************************
[dsoilL,disV,disS,dsoilR,BM,SF] = postprocessor(dis',Eb,Ib,nesL,neB,nesR,leb,L);
d_ = [dsoilL(1:end-1);disV;dsoilR(2:end)]; 
BM(1) = 0; 
BM(end) = 0;
%**************************************************************************
[Fi_z] = phiz(dsoilL,disV,dsoilR,b,H,SpE,nu_,leL,leb,leR,nn1,nn2,nn3,CorX,CorY,CorXn,CorYn);
dFi_z = diff(Fi_z);
dFi_z_dz = dFi_z/dz;
end
return