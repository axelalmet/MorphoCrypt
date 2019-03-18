function Residuals = SelfContactInContactRegionBCs(Ml, Mr, modelParams)
% Boundary conditions specified for self-contact along a region, before the
% region of contact.

% Obtain the relevant model parameters
w0 = modelParams.w0;

Residuals = [Ml(1), ... % S(0) = 0
             Ml(2), ... % x(0) = 0
             Ml(5), ... % G(0) = 0
             Ml(6), ... % theta(0) = 0
             Mr(8) - Ml(8),...  % s(sc) = sc
             Mr(2) - (w0), ... x(sc) = w0
             Mr(6) + 0.5*pi, ... % theta(sc) = -pi/2
             Mr(7)]; % m(sc) = 0

         
end