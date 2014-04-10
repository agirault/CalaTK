function [f g] = functiontest(x)

    f = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
    %f = x(1)^2 + x(2)^2;
end