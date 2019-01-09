function [G0_U,e_new] = G0_Update(ev,e,sm0)
% % % % % % % % % % % % % % % % % % % %
Pa = 100.0;  % Reference pressure (kPa)
Cg = 650;
eg = 2.17;
ng = 0.45;
% 
de = (1 + e).*ev;
e = e - de;
G0_U = Pa * Cg * (eg - e).^2./(1+e).*(sm0/Pa).^ng;
e_new = e;
return