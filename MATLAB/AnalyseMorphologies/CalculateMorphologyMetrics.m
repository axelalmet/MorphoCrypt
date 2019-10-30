function [P, D, H, W, A] = CalculateMorphologyMetrics(rodSol, solMesh)
% Function to calculate various measures of a solution profile of a buckled
% rod.
%
% Input:
% @solution - solution structure, as determined by bvp4c
% @mesh - specified mesh to account for irregular spacings of solution.
%
% Output: 
% P - height of peaks 
% D - depth of invagination
% H - total vertical height of crypt
% W - width of invagination
% A - aspect ratio of invagniation = W/D
% 
% These solutions assume that the rod is invaginating upwards, but you can
% easily adapt the code if it's going downwards. We're also assuming
% half-interval deformations

% Get the relevant solution components
xSol = rodSol.y(2,:);
x = interp1(rodSol.x, xSol, solMesh);

ySol = rodSol.y(3,:);
y = interp1(rodSol.x, ySol, solMesh);

thetaSol = rodSol.y(6,:);
theta = interp1(rodSol.x, thetaSol, solMesh);

D = max(y);
P = min(y);
H = D - P;

% To work out the width is a bit tricky, but we're effectively looking for
% when the region of the rod before theta dips to around -0.5*pi, which is
% based on the self-contact criterion.
index = find(theta < -0.5*pi, 1, 'last');
W = 2*max(x(1:index));

A = W/D;