function PStar = Get2ndCriticalStressValue(parameters, searchInterval)
% Function to determine the critical stress at which a growing,
% inextensible rod buckles.

% Get the first critical growth value
PMin = Get1stCriticalStressValue(parameters);
                        
PCondition = @(p) GetFunctionToDetermine2ndCriticalStress(p, parameters);
                                
PRoot = fzero(PCondition, searchInterval);

if(real(PRoot) ~= PMin)
    PStar = real(PRoot);
else
    PStar = NaN;
end