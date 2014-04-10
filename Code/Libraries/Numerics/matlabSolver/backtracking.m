function [histout,costdata, x] = backtracking(x0,f,g)

%% INITIALISATION
% linesearch parameters (to give in input?)
InitStepSize = 1.0e-03;
MinAllowedStep = 1.0e-04;
MaxItImDecr = 2;
MaxItNonImDecr = 2;
SizeUpFactor = 2.0;
SizeDownFactor = 0.5;
maxit = 100;
backtrackMaxIt = 10;
c1 = 1.d-4;

NbItImDecr    = 0; % Number of iterations with immediate decrease
NbItNonImDecr = 0; % NUmber of iterations without immediate decrease

% current values for x, f(x)
x_cur = x0;
f_cur = feval(f, x_cur); numf = 1; numg = 0;

DesiredStepSize = InitStepSize;

% Hist
it_count = 0;
ithist = zeros(maxit,3);
ithist(1,1) = f_cur;
ithist(1,2) = 0; 
ithist(1,3) = it_count;


%% ITERATIONS
while(it_count <= maxit)
    it_count = it_count+1;
    
    g_cur = feval(g,x_cur); numg = numg + 1;
    Alpha = DesiredStepSize;
    
    %Backtracking loop
    it_backtracking = 0;
    SufficientDecrease = 0;
    while ((Alpha >= MinAllowedStep) & (it_backtracking <= backtrackMaxIt))
        it_backtracking = it_backtracking + 1;
        
        % compute values for x, f(x)
        x_test = x_cur-Alpha*g_cur;
        f_test = feval(f, x_test); numf = numf + 1;
        
        % test armijo condition
        if f_test < f_cur - c1*Alpha*g_cur'*g_cur
            SufficientDecrease = 1;
            x_cur = x_test;
            f_cur = f_test;
            break
        else
            Alpha = Alpha * SizeDownFactor;
        end
    end
           
    if SufficientDecrease == 1        
        if Alpha == DesiredStepSize
            NbItImDecr = NbItImDecr + 1;
            NbItNonImDecr = 0;
        else
            NbItImDecr = 0;
            NbItNonImDecr = NbItNonImDecr + 1;
        end
        
        if NbItImDecr >= MaxItImDecr
            DesiredStepSize = DesiredStepSize * SizeUpFactor;
        end
        if NbItNonImDecr >= MaxItNonImDecr
            DesiredStepSize = DesiredStepSize * SizeDownFactor;
        end
    else
        if Alpha <= MinAllowedStep
            break
        else
            NbItImDecr = 0;
            NbItNonImDecr = NbItNonImDecr + 1;
            DesiredStepSize = Alpha;
        end
    end
    
    % Update Hist
	ithist(it_count,1) = f_cur;
    ithist(it_count,2) = it_backtracking;
	ithist(it_count,3) = it_count; 
    
end

%% OUTPUTS
x = x_cur; 
histout = ithist(1:it_count,:);
costdata = [numf, numg];
