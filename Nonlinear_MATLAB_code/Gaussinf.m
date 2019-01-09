%**************************************************************************
% Gauss Quadrature Method
%**************************************************************************
function [In] = Gaussinf(f)
% 
% 3 - points Gauss Quadrature Method
W1 = 5/9; x1 = sqrt(0.6);
W2 = 8/9; x2 = 0.0;
W3 = 5/9; x3 = -sqrt(0.6);
In = W1*f(x1) + W2*f(x2) + W3*f(x3);
% 
% 4 - points Gauss Quadrature Method
% W1 = 1/2-1/(6*sqrt(1.2)); x1 = sqrt((3+2*sqrt(1.2))/7);
% W2 = 1/2-1/(6*sqrt(1.2)); x2 = -sqrt((3+2*sqrt(1.2))/7);
% W3 = 1/2+1/(6*sqrt(1.2)); x3 = sqrt((3-2*sqrt(1.2))/7);
% W4 = 1/2+1/(6*sqrt(1.2)); x4 = -sqrt((3-2*sqrt(1.2))/7);
% In = W1*f(x1) + W2*f(x2) + W3*f(x3) + W4*f(x4);
% 
return