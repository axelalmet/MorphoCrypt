function [LowerBound, UpperBound] = GetIntervalToFind2ndCriticalStress(stressValues, parameters, tol)
% Function to determine the interval on which we should search for the
% critical growth value that gives rise to buckling in the model of a morphoelastic
% rod on an elastic foundation.

% Get the function we solve to find the critical gamma value.
CriticalStressFunction = @(g) GetFunctionToDetermine2ndCriticalStress(g, parameters); 

% The upper bound should be the first point where the function returns a
% postiive value
upperBoundIndex = find(CriticalStressFunction(stressValues) > tol, 1, 'first');
UpperBound = stressValues(upperBoundIndex);

% The lower bound should be the last point where the function returns a
% negative value. 
lowerBoundIndex = find(CriticalStressFunction(stressValues(1:upperBoundIndex)) < -tol, 1, 'last');
LowerBound = stressValues(lowerBoundIndex);
