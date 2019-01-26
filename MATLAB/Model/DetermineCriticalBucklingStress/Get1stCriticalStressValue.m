function n3sOne = Get1stCriticalStressValue(parameters)
% Function to calculate the 1st critical stress value that gives rise to
% the rod buckling on an infinite domain.

% @param parameters, relevant model parameters (changes depending on
% whether we are considering dimensional or nondimensional model)
% 
% @ returns GammaOne, the calculated gamma function.
    
K = parameters.K;
% Extract the relevant parameter    
n3sOne = 2*sqrt(K);

