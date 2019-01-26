function ExploreSloughingMechanisms
% Script to explore how the sloughing region changes under different
% 'sloughing rate' mechanisms

close all

%% Parameter set-up
L0 = 0.5;
L = L0;
g = 1;
dt = 0.05;
T = 2.5;
tMax = 7;
AStar = 4.0;
sigma = 10;
times = 0:dt:tMax;


%% Homogeneous growth, ramped up sloughing up to growth rate after time t = T
LOld = L;
newLengths = zeros(1, length(times));
newLengths(1) = LOld;

% Iterate using derived numerical scheme for derivatives
for i = 1:(length(times) - 1)
    
    t = times(i);
    
    LNew = LOld + g*dt/(1 + g*t)*(- 0.25*(1 + tanh(sigma*(t - T))) + (L0 - LOld));
    
    LOld = LNew;
    newLengths(i + 1) = LOld;
   
end

% Use 'simplified method', i.e. straight integration with homogeneous parameters.
newLengthsV2 = L0 - cumtrapz(times, 0.25*g*(1 + tanh(sigma*(times - T))))./(1 + g.*times);

% Make sure the numerical scheme isn't too inaccurate.
figure(1)
hold on
plot(times, newLengths)
plot(times, newLengthsV2)


%% Compare Age and Time for homogeneous growth

Ages = times.*(2 + g*times)./(2*(1 + g*times));
figure(2)
hold on
plot(times, Ages)
plot(times, times)
plot(times, 0.5*times)

%% Homogeneous growth, ramped up sloughing up to growth rate after age A = A*
LOld = L;
newLengths = zeros(1, length(times));
newLengths(1) = LOld;
AOld = 0;
gammaOld = 1;

totalLengths = ones(1, length(times));

% Iterate using derived numerical scheme for derivatives
for i = 2:(length(times))
        
    gammaNew = gammaOld + g*dt;
    ANew = 2*AOld + dt - gammaNew/gammaOld*AOld;
    LNew = LOld + g*dt/(gammaOld)*(-0.25*(1 + tanh(sigma*(times(i) - T))) + (L0 - LOld));
    
    gammaOld = gammaNew;
    AOld = ANew;
    LOld = LNew;
    newLengths(i) = LOld;
    totalLengths(i) = LNew*gammaNew;
   
end

figure(1)
hold on
plot(times, newLengths)

figure(2)
hold on
plot(times, totalLengths)
%% Chemical-based heterogeneous growth, age-dependent sloughing

LOld = L;
newLengths = zeros(1, length(times));
newLengths(1) = LOld;
AOld = 0;
gammaOld = 1;
sigma = 0.16;

W = @(S, width) exp(-((S - 0.5*S(end))./width).^2);

S = 0:1e-2:1;
currentS = S;
% Iterate using derived numerical scheme for derivatives
for i = 2:(length(times))
        
    gammaNew = gammaOld.*(1 + dt*W(currentS, sigma));
    
    ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
    
    index = find(S == LOld, 1);
    if (index < length(S))
        LNew = LOld + (-1e-4*(1 + tanh(sigma*(max(AOld) - AStar))));
    else
        LNew = LOld + (-1e-4*(1 + tanh(sigma*(max(AOld) - AStar))));
    end
        
    currentS = cumtrapz(S, gammaNew);
    
    gammaOld = gammaNew;
    AOld = ANew;
    LOld = LNew;
    newLengths(i) = LOld;
   
end

figure(1)
hold on
plot(times, newLengths)
