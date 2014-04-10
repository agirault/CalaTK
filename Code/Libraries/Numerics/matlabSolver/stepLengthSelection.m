function [alpha]=stepLengthSelection(phi_0, phiP_0, alpha_cur, phi_alpha, alpha_prev, phi_prev)
%
% Cubic/quadratic polynomial step length selection
%
% Finds minimizer alpha of the quadratic (first stepsize reduction)
% or the cubic (higher stepsize reduction) polynomial Phi on the interval
% [blow * alpha_cur, bhigh * alpha_cur]
% 

%% INIT
bhigh=.5; blow=.1;
alpha_min = alpha_cur*blow;
alpha_max = alpha_cur*bhigh; 

%% STEP LENTH SELECTION
if nargin == 4 % quadratic model
    %alpha = - phiP_0/(2 * alpha_cur*(phi_alpha - phi_0 - phiP_0) );
    alpha = - phiP_0*alpha_cur^2/(2*(phi_alpha - phi_0 - phiP_0*alpha_cur) );
    if alpha < alpha_min alpha = alpha_min; end
    if alpha > alpha_max alpha = alpha_max; end
    
else % cubic model
    %a=[alpha_cur^2, alpha_cur^3; alpha_prev^2, alpha_prev^3];
    %b=[phi_alpha; phi_prev]-[phi_0 + phiP_0*alpha_cur; phi_0 + phiP_0*alpha_prev];
    %c=a\b;
    %alpha=(-c(1)+sqrt(c(1)*c(1) - 3 *c(2) *phiP_0))/(3*c(2));
    
    M1 = [alpha_prev^2, -alpha_cur^2; -alpha_prev^3, alpha_cur^3];
    M2 = [phi_alpha - phi_0 - phiP_0*alpha_cur; phi_prev - phi_0 - phiP_0*alpha_prev];
    M = M1 * M2 / (alpha_prev^2 * alpha_cur^2 * (alpha_cur - alpha_prev));
    a = M(1);
    b = M(2);
    alpha = (-b + sqrt(b*b - 3*a*phiP_0))/(3*a);
    
    if alpha < alpha_min alpha = alpha_min; end
    if alpha > alpha_max alpha = alpha_max; end
end