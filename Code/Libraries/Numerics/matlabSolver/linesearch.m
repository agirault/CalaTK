function [histout,costdata, x] = linesearch(x0,f,g,tolPos,tolGrad,maxit)
%
% LINESEARCH finds the minimum of an objective function using
%	- Steepest descent Algorithm
%	- Armijo as a condition for sufficient decrease (could try Wolf?)
%	- Step length selection
% 
%   INPUT:
%       - x0 = initial iterate position
%       - f = objective function,
%       - g = grad f (COLUMN vector)
%       - tolPos = termination criterion norm(dx) < tolPos
%       (optional) default = 1.d-6
%       - tolGrad = termination criterion norm(grad) < tolGrad
%       (optional) default = 1.d-6
%       - maxit = maximum iterations
%       (optional) default = 1000
%
%   OUTPUT:
%       - x = solution (minimum found)
%       - histout = iteration history. Each row of histout is :
%       [norm(grad), f, number of step length reductions, iteration count]
%       - costdata = [num f, num grad]
%
% Requires: stepLengthSelection.m
%

%% INITIALISATION
% linesearch parameters
if nargin < 6
    maxit = 1000; 
end
if nargin < 5
    tolGrad = 1.d-6;
end
if nargin < 4
    tolPos = 1.d-6;
end

c1 = 1.d-4;

% current values for x, f(x) and g(x)
dx = 1.d10;
x_cur = x0; 
f_cur = feval(f,x_cur); numf = 1;
g_cur = feval(g,x_cur); numg = 1;
%[f_cur,g_cur] = feval(f,x_cur); numf = 1; numg = 1;

% Hist
it_count = 0;
ithist = zeros(maxit,4);
ithist(1,1) = norm(g_cur);
ithist(1,2) = f_cur;
ithist(1,3) = 0; 
ithist(1,4) = it_count;

%% ITERATIONS
while((norm(g_cur) > tolGrad) & (dx > tolPos) & (it_count <= maxit))
    it_count = it_count+1;
    
    % Compute alpha0 and the values to test Armijo condition
	alpha_cur   =   min(1,100/(1+norm(g_cur)));   %fixup for very long steps
    x_test      =	x_cur-alpha_cur*g_cur;
    phi_alpha   =	feval(f,x_test); numf=numf+1;
    phi_0       =	f_cur;
    phiP_0      =	-g_cur'*g_cur;
    
    % Test Armijo Condition
    it_armijo = 0;
	while( phi_alpha > (phi_0 + c1*alpha_cur*phiP_0) )
        
        % test # of iterations
        it_armijo = it_armijo+1;
		if(it_armijo > 10) 
            disp(' Armijo error in steepest descent')
            histout = ithist(1:it_count,:); costdata=[numf, numg, numh];
            return;
        end
        
        % compute new alpha with step length selection
        if it_armijo == 1 % quadratic model
            alpha = stepLengthSelection(phi_0, phiP_0, alpha_cur, phi_alpha);
        else % cubic model
            alpha = stepLengthSelection(phi_0, phiP_0, alpha_cur, phi_alpha, alpha_prev, phi_prev);
        end
         
        % store previous values
        phi_prev    =   phi_alpha;
        alpha_prev  =   alpha_cur;
        alpha_cur   =   alpha;
        
        % calculate new value to test for Armijo
		x_test      =   x_cur-alpha_cur*g_cur;
		phi_alpha   =   feval(f,x_test); numf = numf+1;

    end
    
    % Update new values x, f(x) and g(x)
    dx = norm( x_test - x_cur );
    x_cur = x_test;
    f_cur = phi_alpha;
    g_cur = feval(g,x_cur); numg = numg+1;
    %[f_cur,g_cur] = feval(f,x_cur); numf = numf+1; numg = numg+1;
    
    % Update Hist
	ithist(it_count,1) = norm(g_cur);
    ithist(it_count,2) = f_cur;
    ithist(it_count,3) = it_armijo;
	ithist(it_count,4) = it_count; 
end

%% OUTPUTS
x = x_cur; 
histout = ithist(1:it_count,:);
costdata = [numf, numg];
