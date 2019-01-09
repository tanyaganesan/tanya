%**************************************************************************
% PROGRAM FOR BEAME ON NON-LINEAR FOUNDATION USING MODIFIED VLASOV'S MODEL
%**************************************************************************
% function Main_Prog(N_,L,M_ratio,COV_Es,Corr_Dist_Es,f1,f2,f3,f4,f5,f6)
clc                         % Clears the command windows
clear all                   % Clear all stored variables
Tf = datestr(now,13)        % Prints current computer time
%**************************************************************************
% Deterministic Parameters for Beam & soil
%**************************************************************************
L = 4;                      % Length of the Beam (m)
b = 0.5;                    % Width of the Beam (m)
d = 0.5;                    % Depth of the Beam (m)
Ab = d * b;                 % C/s area of beam (m^2)
Eb = 25e6;               % Elastic Modulus of the Beam (kPa)
Ib = (1/12)*b*d^3;          % Moment of Inertia of Beam (m^4)
H = 20;                     % Depth of Soil (m)
nu_ = 0.3;                 % Poission's ratio of soil
Es = Eb/100;
%**************************************************************************
% Discretization of Soil Domain
%**************************************************************************
nesL = 50;              % Number of elements : soil (left) 
neB  = 50;              % Number of elements : beam 
nesR = 50;              % Number of elements : soil (right) 
n    = 80;              % Division in Y direction 
m = nesL + neB + nesR;  % Division in X direction
%**************************************************************************
% Size of Soil Domain
%**************************************************************************
Lx_l = 5*L;             % Lateral extent of soil present at left side of beam (m)
Lx_r = 5*L;             % Lateral extent of soil present at right side of beam (m)
Lx  = Lx_l + L + Lx_r;  % Horizontal dimension of the soil domain (m)
Ly  = H;                % Depth of soil (m)
dx  = L/neB;            % Beam element length (m)
dxl = Lx_l/nesL;        % Soil element length : soil (left) (m)
dxr = Lx_r/nesR;        % Soil element length : soil (right) (m)
dz  = H/n;              % Soil element length along depth (m)
%**************************************************************************
[CODRX,CORDY,CorXn,CorYn] = coord(dxl,dx,dxr,dz,H,Lx_l,L,Lx_r,Ly);
% [SpE] = RandomField(m,n,Es,0.5,2,L,Ly,CorX,CorY);                  
gamma_z = 18;                                          % Unit weight of soil (kN/m^3)
Pa = 100;                                              % Reference stress (kPa)
e0 = 0.2;                                              % Initial void ratio of soil
sz0 = CORDY*gamma_z + 2*gamma_z;                       % Geostatic sress [surcharge = 2 x gamma_z] (kPa)
sm0 = 1/3*(sz0 + 0.5*sz0 + 0.5*sz0);                   % Mean effective stress (kPa)
G0 = 650 * (2.17 - e0)^2/(1+e0)*Pa^(1-0.45)*(sm0).^0.45;           % Initial shear modulus of soil : sand(kPa)
% G0 = 650 * (2.17 - e0)^2/(1+e0)*(sm0).^0.45*OCR^0.5; % Initial shear modulus of soil : clay(kPa)
SpE = 2*G0*(1+nu_);                                    % Initial Young's modulus of soil (kPa)
e = e0*ones(n,m);
%**************************************************************************
% Loads on Beam
%**************************************************************************
Q =  -50;  % Value of uniformly distributed Load (kN/m) | -ve sign : downward
V  = -00;  % Value of concentrated load (kN) | -ve sign : downward
Mo = -00;  % Value of concentrated moment (kNm) | +ve sign : Anticlockwise
l_pos = [-1;0.5*L;-1]; % Position of concentrated loads
%**************************************************************************
% gm_oct_prv = 1; gm_oct_new = 2; JJ = 0; gm_oct(1,nesL) = 1;
for JJ = 1:10   % <------ Iteration for non-linearity starts here
% while abs(gm_oct_new - gm_oct_prv) > 1e-6
JJ
% JJ = JJ + 1;    
%**************************************************************************
[d_,disV,BM,SF,nn2,k_,t_,z,Fz,dFzdz] = analysis(nesL,neB,nesR,m,n, ...
           SpE,L,Lx,Ly,nu_,H,b,Eb,Ib,Q,V,Mo,l_pos,dxl,dx,dxr,dz,Lx_l,Lx_r);
%**************************************************************************
% d2v_dx2 = diff(diff(disV)/dx)/dx;
cont_pr = -(b/b)*(k_(:,end).*disV - t_(:,end).*BM/(Eb*Ib));
%**************************************************************************
% gm_oct_prv = gm_oct(1,nesL);
[gm_oct,tao_xz,sig_mn,ezz,exz,ev] = strss_strain(SpE,CORDY,Fz(:,end), ...
            dFzdz(:,end),n,m,nesL,neB,nesR,d_,disV,dxl,dx,dxr,nu_,gamma_z);
% gm_oct_new = gm_oct(1,nesL);
gm_oct_(JJ,1) = gm_oct(1,nesL);
for I1 = 1:n
    for J1 = 1:m
        [G_G0_] = mod_reduction(0.001,gm_oct(I1,J1));
        G_G0(I1,J1) =   G_G0_;
    end
end
[G0_U, e_new] = G0_Update(ev,e,sig_mn);
SpE = 2* G_G0.* G0_U* (1 + nu_);
e = e_new;
%**************************************************************************
end % <------ Iteration for non-linearity ends here
%**************************************************************************
gm_oct_

% [G0_U] = G0_Update(SpE,nu_,ezz,exz,ev)
figure
post_processor_plot(CODRX,CORDY,gm_oct);

%break

X0 = linspace(0,L,nn2);
subplot 411
hold on
plot(X0,disV,'-r','linewidth',2)
% axis([0 L -max(abs(disV))-0.01 0])
xlabel('x (m)')
ylabel('deflection (m)')
grid on
subplot 412
hold on
plot(X0,BM ,'-r','linewidth',2)
% axis([0 L -max(abs(BM))-2 0])
xlabel('x (m)')
ylabel('bending moment (kN-m)')
grid on
subplot 413
hold on
plot(X0,SF,'-r','linewidth',2)
% axis([0 L -max(abs(SF))-2 max(abs(SF))+2])
xlabel('x (m)')
ylabel('shear force (kN)')
grid on
subplot 414
hold on
plot(X0,cont_pr,'-r','linewidth',2)
% axis([0 L 0 max(abs(cont_pr))+2])
xlabel('x (m)')
ylabel('contact pr (kPa)')
grid on
Tf2 = datestr(now,13)                       % Prints current computer date & time