close all;
clear all;

x0 = [5000;5000]
f = @functiontest;
g = @gradienttest;
tolPos  = 1.d-4;
tolGrad = 1.d-2;


[histout,costdata,x] = linesearch( x0, f, g, tolPos, tolGrad)
%[histout,costdata,x] = backtracking( x0, f, g)

%       - x = solution (minimum found)
%       - histout = iteration history. Each row of histout is :
%       [norm(grad), f, number of step length reductions, iteration count]
%       - costdata = [num f, num grad]