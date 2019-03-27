function ExploreSloughingMechanisms
% Script to explore how the sloughing region changes under different
% 'sloughing rate' mechanisms

clear all
%% Parameter set-up
L0 = 1;
L = L0;
g = 1;
dt = 0.05;
T = 5;
tMax = 20.0;
AStar = 5;
sigma = 10;
times = 0:dt:tMax;

%% Homogeneous, time-dependent sloughing

LOld = 1;
newLengths = zeros(1, length(times));
newLengths(1) = 1;
AOld = 0;
gammaOld = 1;
g = 1;
sigma = 0.16;

W = @(S, width) exp(-((S)./width).^2);

S = 0:1e-3:1;
currentS = S;

options = optimset('Display','off');

% Iterate using derived numerical scheme for derivatives
for i = 2:(length(times))
    
    gammaNew = gammaOld + g*dt;
    
    ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
    
    t = times(i);
    
    sloughedAmount = trapz(0:1e-2:t, 0.5*(1 + tanh(10*((0:1e-2:t) - T))) );
    
    currentArcLength = gammaNew.*S;
    
    if (sloughedAmount > 1e-4)
        % Define new functionx
        sloughLfunct = @(l) sloughedAmount - (currentArcLength(end) - interp1(S, currentArcLength, l));
        
        sloughL = min(fsolve(sloughLfunct, [0 1], options));
        newLengths(i) = sloughL;
        
    else
        newLengths(i) = 1;
    end
    
    gammaOld = gammaNew;
    AOld = ANew;
    
end

figure
plot(times, newLengths)

%% Chemical signalling, time-dependent sloughing
LOld = 1;
newLengths = zeros(1, length(times));
newLengths(1) = 1;
AOld = 0;
gammaOld = 1;
sigma = 0.16;

W = @(S, width) exp(-((S)./width).^2);

S = 0:1e-3:1;
currentS = S;

options = optimset('Display','off');

% Iterate using derived numerical scheme for derivatives
for i = 2:(length(times))
    
    gammaNew = gammaOld.*(1 + dt*W(currentS, sigma));
    
    ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
    
    t = times(i);
    
    netGrowth = trapz(S, gammaOld.*W(currentS, sigma));
    
    sloughedAmount = trapz(0:1e-2:t, 0.5*netGrowth*(1 + tanh((10*((0:1e-2:t) - T)))));
    
    currentS = cumtrapz(S, gammaNew);
    
    if (sloughedAmount > 1e-4)
        % Define new functionx
        sloughLfunct = @(l) sloughedAmount - (currentS(end) - interp1(S, currentS, l));
        
        sloughL = min(fsolve(sloughLfunct, [0 1], options));
        newLengths(i) = sloughL;
        
    else
        newLengths(i) = 1;
    end
    
    gammaOld = gammaNew;
    AOld = ANew;
    
end

%%
figure
plot(times, newLengths, '-^')

%% Chemical signalling, age-dependent sloughing
LOld = 1;
newLengths = zeros(1, length(times));
newLengths(1) = LOld;
AOld = 0;
gammaOld = 1;
sigma = 0.16;

W = @(S, width) exp(-((S)./width).^2);

S = 0:1e-3:LOld;
currentS = S;

sloughedAmount = zeros(1, length(times));

options = optimset('Display','off');

% Iterate using derived numerical scheme for derivatives
for i = 2:(length(times))
    
    gammaNew = gammaOld.*(1 + dt*W(currentS, sigma));
    
    ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
    
    sloughedAmount(i) = trapz(S, ones(1, length(S)).*(ANew > AStar).*(AOld.*W(currentS, sigma) < 0.01));
    
    gammaOld = gammaNew;
    AOld = ANew;
    
    if (sloughedAmount(i) > 1e-4)
        
        break
        
    else
        newLengths(i) = 1;
    end
    
    figure(1)
    hold on
    plot(S, ANew)
    
end

%%

% Get the new sloughing boundary

% Define new functionx
sloughLfunct = @(l) sloughedAmount(i) - (currentS(end) - interp1(S, currentS, l));

sloughL = fsolve(sloughLfunct, 1, options);
newLengths(i) = sloughL;

L = sloughL;

SOld = S;
gammaPrev = gammaOld;
APrev = AOld;

S = [SOld(1:(end - 1))./(SOld(end - 1))*L, L, 1];

gammaOld = [gammaPrev(1:(end - 1)), gammaPrev((end - 1):end)];
currentS = cumtrapz(S, gammaOld);

AOld = [APrev(1:(end - 1)), APrev((end - 1):end)];

%%

% Iterate now with the evolving sloughing boundary
for j = (i + 1):length(times)
    gammaNew = gammaOld.*(1 + dt*W(currentS, sigma));
    
    ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
    
    sloughedAmount(j) = trapz(S, ones(1, length(S)).*(ANew > AStar).*(AOld.*W(currentS, sigma) < 0.01));
    
    LOld = newLengths(j - 1);
    LIndex = find(S == LOld, 1);
    
    % Define function to determine new boundary
    sloughLfunct = @(l) sloughedAmount(j) - (currentS(end) ...
        - interp1(S([1:LIndex, (LIndex + 2):end]), ...
        currentS([1:LIndex, (LIndex + 2):end]), l));
    
    sloughL = fsolve(sloughLfunct, newLengths(j - 1), options);
    newLengths(j) = sloughL;
    
    S(LIndex:(LIndex + 1)) = sloughL;
    
    gammaOld = gammaNew;
    AOld = ANew;
    
    currentS = cumtrapz(S, gammaOld);
    
    figure(1)
    hold on
    plot(S, ANew)
    
end

figure(2)
hold on
plot(times, newLengths)
plot(times, sloughedAmount)
