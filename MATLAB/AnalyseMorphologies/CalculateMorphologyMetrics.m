function [P, I, H, W, A] = CalculateMorphologyMetrics(solution, mesh)
% Function to calculate various measures of a solution profile of a buckled
% rod.
%
% Input:
% @solution - solution structure, as determined by bvp4c
% @mesh - specified mesh to account for irregular spacings of solution.
%
% Output: 
% P - height of peaks 
% I - depth of invagination
% H - total vertical height of crypt
% W - width of invagination
% A - aspect ratio of invagnation = W/I