function Residuals = SelfContactOutContactRegionBCs(Ml, Mr, modelParams)
% Boundary conditions specified for self-contact along a region, before the
% region of contact.

% Obtain the relevant model parameters
w0 = modelParams.w0;
y0 = modelParams.y0;
L = modelParams.L;
% sc = modelParams.sc;

Residuals = [Ml(2) - (w0), ... x(sc) = w0
             Ml(6) + 0.5*pi, ... % theta(sc) = -pi/2
             Ml(7), ... % m(sc) = 0
             Mr(1) - L,...  % s(L) = L
             Mr(2) - 0.5, ... x(L) = 0.5
             Mr(3) - y0, ... y(L) = y0
             Mr(6)]; % theta(L) = 0

         
end