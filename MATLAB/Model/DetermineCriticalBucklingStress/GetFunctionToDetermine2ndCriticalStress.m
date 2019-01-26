function StressFunction = GetFunctionToDetermine2ndCriticalStress(P, parameters)
% Function of which the zeros correspond to the critical growth values for
% which buckling occurs.
%
% @param gamma, the input value
% @param parameters, the relevant model parameters
% @param isDimensional, 1 if using for dimensional model, 0 otherwise
% (determines how we define function). 
%
% Return values of function at gamma.
%
% Author: Axel Almet
% Created on: 17/04/17
% Last modified: 17/04/17

K = parameters.K;
    
% Define the ODE coefficients
a = @(p) p./2;
b = @(p) sqrt(K);


w1 = @(p) sqrt(a(p) + sqrt((a(p)).^2 - (b(p)).^2));
w2 = @(p) sqrt(a(p) - sqrt((a(p)).^2 - (b(p)).^2));
                    
StressFunction = a(P).*sin(w1(P)).*sin(w2(P)) + b(P).*cos(w1(P)).*cos(w2(P)) - b(P);