%**************************************************************************
% Determines beam displacement, bending moment & shear force
%**************************************************************************
function [d_soil_L,dis_b_V,dis_b_S,d_soil_R,bm,sf] = postprocessor(d,Eb,Ib,nesL,neB,nesR,leb,L)
nn = neB + 1;
element2 = [1:1:(nn-1); 2:1:nn]'; % Nodal Connectivity of beam
% % 
d_soil_L = d(1:nesL+1);
d_beam = d(nesL+1:(nesL+1)+2*(neB+1)-1);
d_soil_R = [d((nesL+1)+2*(neB+1)-2);d((nesL+1)+2*(neB+1):end)];
NQ  = max(size(d_beam));
for I = 1:NQ/2
disV(I) = d_beam(2*I-1);
disS(I) = d_beam(2*I);
end
% % 
dis_b_V = disV';
dis_b_S = disS';
% Calculation of Bending Moment
BM = zeros(nn,1);
SF = zeros(nn,1);
x1 = -sqrt(0.6);
x2 = sqrt(0.6);
[Nb1,d2Nb1,d3Nb1,Ns1] = shape_fn(leb,x1);
BM1 = d_beam([1:2:(2*nn-3); 2:2:2*nn-2 ; 3:2:2*nn-1 ; 4:2:2*nn]')*0.5*d2Nb1';
SF1 = d_beam([1:2:(2*nn-3); 2:2:2*nn-2 ; 3:2:2*nn-1 ; 4:2:2*nn]')*0.5*d3Nb1';
[Nb2,d2Nb2,d3Nb2,Ns2] = shape_fn(leb,x2);
BM2 = d_beam([1:2:(2*nn-3); 2:2:2*nn-2 ; 3:2:2*nn-1 ; 4:2:2*nn]')*0.5*d2Nb2';
SF2 = d_beam([1:2:(2*nn-3); 2:2:2*nn-2 ; 3:2:2*nn-1 ; 4:2:2*nn]')*0.5*d3Nb2';
M = Eb * Ib * [BM1'; BM2'];
S = Eb * Ib * [SF1'; SF2'];
% %  
for e2 = [1:neB];  
    ss2 = element2(e2,:);
    BM(ss2,1) = BM(ss2,1) + M(:,e2);
    SF(ss2,1) = SF(ss2,1) + S(:,e2);
end
% %  
X0 = linspace(0,L,nn);
bm = interp1(X0(2:end-1)',BM(2:end-1),X0','spline','extrap');
sf = interp1(X0(2:end-1)',SF(2:end-1),X0','spline','extrap');
return