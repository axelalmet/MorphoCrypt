function InvestigateGrowthLawsWithMechanicalFeedback
% Script where we can test different growth hypotheses.
close all

tMax = 10;
dt = 0.05;

times = 0:dt:tMax;

W = @(S, width) exp(-((S)./width).^2);

initArcLength = 0:5*1e-3:1;
currentArcLength = initArcLength;
initGamma = ones(1, length(initArcLength));

sigma = 2*10/125;

gammaSols = cell(length(times), 1);
gammaSols{1} = [initArcLength; initGamma];

mu = 10;
Es = 1;
ns = -0.4;

outputDirectory = '../Solutions/UnbuckledRod/';
outputValues = 'linearmechanofeedback_Es_1_K_infty_ns_-0p6_mu_1';

%%

mu = 1;

outputValues = 'linearmechanofeedback_Es_1_K_infty_ns_-0p4_mu_1';

gammaOld = initGamma;

for i = 2:length(times)
    
    gammaNew = gammaOld.*(1 + dt*(W(initArcLength, sigma) + mu.*(Es.*((1 - gammaOld)./(gammaOld)) - ns)));
    
    gammaSols{i} = [initArcLength; gammaNew; W(initArcLength, sigma) + mu.*(Es.*((1 - gammaOld)./(gammaOld)) - ns); Es.*(1 - gammaNew)./gammaNew];
    
    %     figure(1)
    %     hold on
    %     plot(initArcLength, gammaNew)
    %
    %     figure(2)
    %     hold on
    %     plot(initArcLength, (gammaNew - gammaOld)./(dt*gammaOld))
    
    gammaOld = gammaNew;
    
    
    
end

save([outputDirectory, 'gammasols_', outputValues, '.mat'], 'gammaSols'); % Solutions
save([outputDirectory, 'gammasoltimes_', outputValues, '.mat'], 'times'); % Solutions
%% Plot gamma for set values of w

mu = 10;

outputValues = 'linearmechanofeedback_Es_1_K_infty_ns_-0p4_mu_10';

wValues = 0:0.1:1;

GammaFixedSSols = cell(length(wValues), 1);

for i = 1:length(wValues)
    w = wValues(i);
    
    gamma = mu*Es/(mu*(Es + ns) - w) + (mu*ns - w)./(mu*(Es + ns) - w).*exp(-(mu*(Es + ns) - w).*times);
    
    gammaIncremental = (w - mu*ns).*(mu*(Es + ns) - w).*exp(-(mu*(Es + ns) - w).*times)./(mu*Es + (mu*ns - w).*exp(-(mu*(Es + ns) - w).*times));
    
    figure(3)
    hold on
    plot(times, gamma)
    
    figure(4)
    hold on
    plot(times, gammaIncremental)
    
    GammaFixedSSols{i} = [gamma; gammaIncremental];
    
end

save([outputDirectory, 'gammasols_fixedwnt_', outputValues, '.mat'], 'GammaFixedSSols'); % Solutions
save([outputDirectory, 'times_fixedwnt_', outputValues, '.mat'], 'times'); % Solutions
save([outputDirectory, 'wntvalues_fixedwnt_', outputValues, '.mat'], 'wValues'); % Solutions



%%
mu = 1;

outputValues = 'localmechanofeedback_Es_1_K_infty_ns_-0p4_mu_1';

gammaOld = initGamma;

for i = 2:length(times)
    
    gammaNew = gammaOld.*(1 + dt*(W(initArcLength, sigma) + mu.*(Es.*(1 - gammaOld)./(gammaOld) - ns).*((Es.*(1 - gammaOld)./(gammaOld)) < ns)));
    
    gammaSols{i} = [initArcLength; gammaNew; (gammaNew - gammaOld)./(dt*gammaOld); Es.*(1 - gammaNew)./gammaNew];
    
        figure(1)
        hold on
        plot(initArcLength, gammaNew)
    
        figure(2)
        hold on
        plot(initArcLength, (gammaNew - gammaOld)./(dt*gammaOld))
    
    gammaOld = gammaNew;
    
    
    
end

save([outputDirectory, 'gammasols_', outputValues, '.mat'], 'gammaSols'); % Solutions
save([outputDirectory, 'gammasoltimes_', outputValues, '.mat'], 'times'); % Solutions

%% Plot gamma for set values of w

mu = 1;

wValues = 0:0.1:1;

GammaFixedSSols = cell(length(wValues), 1);

for i = 1:length(wValues)
    w = wValues(i);
    
    gamma = exp(w*times).*(times < 1/w*log(Es/(Es + ns))) + (times >= 1/w*log(Es/(Es + ns))).*(mu*Es/(mu*(Es + ns) - w) - w*Es./((Es + ns)*(mu*(Es + ns) - w)).*exp(-(mu*(Es + ns) - w).*(times - 1/w*log(Es/(Es + ns)))));
    
    gammaIncremental = w.*(times < 1/w*log(Es/(Es + ns))) + ...
        (times >= 1/w*log(Es/(Es + ns))).*(w*Es).*(mu*(Es + ns) - w).*exp(-(mu*(Es + ns) - w).*(times - 1/w*log(Es/(Es + ns))))./(mu*Es*(Es + ns) - w*Es.*exp(-(mu*(Es + ns) - w).*(times - 1/w*log(Es/(Es + ns)))));
    
    figure(3)
    hold on
    plot(times, gamma)
    
    figure(4)
    hold on
    plot(times, gammaIncremental)
    
    GammaFixedSSols{i} = [gamma; gammaIncremental];
    
end
% 
outputValues = 'localmechanofeedback_Es_1_K_infty_ns_-0p4_mu_1';
% 
save([outputDirectory, 'gammasols_fixedwnt_', outputValues, '.mat'], 'GammaFixedSSols'); % Solutions
save([outputDirectory, 'times_fixedwnt_', outputValues, '.mat'], 'times'); % Solutions
save([outputDirectory, 'wntvalues_fixedwnt_', outputValues, '.mat'], 'wValues'); % Solutions

