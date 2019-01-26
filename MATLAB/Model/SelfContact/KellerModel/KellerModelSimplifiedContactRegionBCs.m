function Residuals = KellerModelSimplifiedContactRegionBCs(Ml, Mr)
% Boundary conditions specified for self-contact along a region.

Residuals = [Ml(1), Ml(2), Ml(3), Mr(2), Ml(5), Ml(6), Mr(6) - 0.5*pi, Mr(7)]; % theta3(L) = 0
         
end