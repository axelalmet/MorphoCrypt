function InvestigateGrowthLaws
% Script where we can test different growth hypotheses.
close all

tMax = 5;
dt = 1;

times = 0:dt:tMax;

W = @(S, width) exp(-((S - 0.5*S(end))./width).^2);

initArcLength = 0:1e-2:1;
currentArcLength = initArcLength;
initGamma = 1;

sigma = 10/125;
%% Let's first test Wnt-only signalling, where Wnt is parametrised by S0

gammaOld = initGamma;

figCount = 1;
% Test linear growth, Wnt parametrised by initial arc length
for i = 1:length(times)
    
    gammaNew = gammaOld + dt*W(initArcLength, sigma);
    
    figure(figCount)
    hold on
    plot(initArcLength, gammaNew)
    
    gammaOld = gammaNew;
end

figure(figCount)
hold on
title('Linearly increasing growth, parametrised by S_0')

gammaOld = initGamma;

figCount = figCount + 1;

% Test exponential growth, Wnt parametrised by initial arc length
for i = 1:length(times)
    
    gammaNew = gammaOld.*(1 + dt*W(initArcLength, sigma));
    
    figure(figCount)
    hold on
    plot(initArcLength, gammaNew)
    
    gammaOld = gammaNew;
end

figure(figCount)
hold on
title('Exponentially increasing growth, parametrised by S_0')

gammaOld = initGamma;
figCount = figCount + 1;
% Test compartment-based growth
for i = 1:length(times)
    
    gammaNew = gammaOld.*(1 + dt.*(initArcLength < (0.5 + sigma).*(initArcLength > 0.5 - sigma)));
    
    figure(figCount)
    hold on
    plot(initArcLength, gammaNew)
    
    gammaOld = gammaNew;
end

figure(figCount)
hold on
title('Compartment-based growth, parametrised by S_0')

%%
gammaOld = initGamma;

% Test linear growth, Wnt parametrised by initial arc length
for i = 1:length(times)
    
    gammaNew = gammaOld + dt*W(currentArcLength, sigma);
    currentArcLength = cumtrapz(initArcLength, gammaNew);
    
    figure(4)
    hold on
    plot(initArcLength, gammaNew)
    
    gammaOld = gammaNew;
end

figure(4)
hold on
title('Linearly increasing growth, parametrised by s')

gammaOld = initGamma;

% Test exponential growth, Wnt parametrised by initial arc length
for i = 1:length(times)
    
    gammaNew = gammaOld.*(1 + dt*W(currentArcLength, sigma));
    currentArcLength = cumtrapz(initArcLength, gammaNew);
    
    figure(5)
    hold on
    plot(initArcLength, gammaNew)
    
    gammaOld = gammaNew;
end

figure(5)
hold on
title('Exponentially increasing growth, parametrised by s')

gammaOld = initGamma;

% Test compartment based growth, parametrised by current arc length
for i = 1:length(times)
    
    gammaNew = gammaOld.*(1 + dt.*(currentArcLength < (0.5*currentArcLength(end) + sigma) ...
                    .*(currentArcLength > 0.5*currentArcLength(end) - sigma)));
    currentArcLength = cumtrapz(initArcLength, gammaNew);
        
    figure(6)
    hold on
    plot(initArcLength, gammaNew)
    
    gammaOld = gammaNew;
end

figure(6)
hold on
title('Compartment-based growth, parametrised by s')

