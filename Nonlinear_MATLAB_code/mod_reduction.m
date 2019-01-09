function [G_Gmax] = mod_reduction(gm_r,gm_)
% Generates modulus reduction curves

f = 0.98;
g = 0.25;
% gm_r = 0.001;
% gm_ = linspace(1e-6,.01,10000)';
% G_Gmax = 1 - f*(tao_s/tao_s_max)^g;
% gm_ = 5.0172e-008
% for I = 1:10000
X = fzero(@(x) x + f*x*x^(g-1)*(gm_/gm_r)^g - 1,1);
% end
% semilogx(gm_,X)
G_Gmax = X;


% Hardin Model
% G_Gmax = 1./(1 + gm_/gm_r);

% Ishibashi and Zhang (1993) model
% Ref: Ishibashi and Zhang (1993)
% PI : Plasticity Index of Soil
% sm : Mean effective stress
% if PI == 0
%     n = 0.00;                   % PI = 0
% else if 0 < PI <= 15
%         n = 3.37e-6*PI^1.404;     % 0 < PI <= 15
%     else if 15 < PI <= 70
%             n = 7.00e-7*PI^1.976;     % 15 < PI <= 70
%         else if PI > 70
%                 n = 2.70e-5*PI^1.115;     % PI > 70
%             end
%         end
%     end
% end
% % 
% K = 0.5*(1 + tanh(log( ((0.000102 + n)./gm_).^0.492)));
% % 
% A = 0.272*(1 - tanh(log( (0.000556./gm_).^0.4)))*exp(-0.0145*PI^1.3);
% % 
% G_Gmax = K .* sm.^A;
return