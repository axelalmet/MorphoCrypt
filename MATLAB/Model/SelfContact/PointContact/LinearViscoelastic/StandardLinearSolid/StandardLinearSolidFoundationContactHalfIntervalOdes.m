function Odes = StandardLinearSolidFoundationContactHalfIntervalOdes(x, M, region, MOld, modelParams)

L = modelParams.L;
gamma = modelParams.gamma;
Eb = modelParams.Eb;
y0 = modelParams.y0;
K = modelParams.K;
beta = modelParams.beta;
L = modelParams.L;
nu = modelParams.nu;
dt = modelParams.dt;

contIndex = find(MOld.x == 1);
xOld = MOld.x([1:contIndex, (contIndex + 2):end]);

SOld = MOld.y(1,:);
XOld = MOld.y(2,:);
YOld = MOld.y(3,:);
POld = modelParams.P;

    function dMdS = ContactEqns(x, M, region)
        
        % Define the state variables
        S = M(1,:);
        X = M(2,:);
        Y = M(3,:);
        F = M(4,:);
        G = M(5,:);
        theta = M(6,:);
        m = M(7,:);
        sc = M(8,:);
        
        % Interpolate the parameters if they are non-constant
        if (length(gamma) > 1)
            gammaOld = gamma([1:contIndex, (contIndex + 2):end]);
            gamma = interp1(xOld, gammaOld, x);
        end
        
        if (length(Eb) > 1)
            EbOld = Eb([1:contIndex, (contIndex + 2):end]);
            Eb = interp1(xOld, EbOld, x);
        end
        
        if (length(POld) > 1)
            POld = POld([1:contIndex, (contIndex + 2):end]);
            POld = interp1(xOld, POld, x);
        end
        
        if (length(SOld) > 1)
            SOld = SOld([1:contIndex, (contIndex + 2):end]);
            SOld = interp1(xOld, SOld, x);
        end
        
        if (length(XOld) > 1)
            XOld = XOld([1:contIndex, (contIndex + 2):end]);
            XOld = interp1(xOld, XOld, x);
        end
        
        if (length(YOld) > 1)
            YOld = YOld([1:contIndex, (contIndex + 2):end]);
            YOld = interp1(xOld, YOld, x);
        end
        
        
        DeltaOld = sqrt((XOld - SOld).^2 + (YOld).^2);
        Delta = sqrt((X - S).^2 + (Y).^2);
        
        % Update the viscoelastic stress in time
        if (nu == 0)
            
            P = (1 - beta)*K.*(Delta - y0);
            
        else
            
            P = (1 - dt*beta/nu)^(-1).*(POld + dt*beta*(1 - beta)*K/nu.*(Delta - y0) ...
                     + (K/nu).*(Delta - DeltaOld));
            
        end
        
        switch region
            case 1 % ODEs for [0, sc1]
                
                l = sc;
                
                
                dSdS = (l).*ones(1, length(S));
                dxdS = (l).*gamma.*cos(theta);
                dydS = (l).*gamma.*sin(theta);
                dFdS = (l).*gamma.*P./Delta.*(X - S);
                dGdS = (l).*gamma.*P./Delta.*(Y);
                dthetadS = (l).*gamma.*(m./Eb);
                dmdS = (l).*gamma.*(F.*sin(theta) - G.*cos(theta));
                dscdS = zeros(1, length(S));
                
            case 2 % ODEs for [sc1, L/2]
                
                l = L - sc;
                
                dSdS = (l).*ones(1, length(S));
                dxdS = (l).*gamma.*cos(theta);
                dydS = (l).*gamma.*sin(theta);
                dFdS = (l).*gamma.*P./Delta.*(X - S);
                dGdS = (l).*gamma.*P./Delta.*(Y);
                dthetadS = (l).*gamma.*(m./Eb);
                dmdS = (l).*gamma.*(F.*sin(theta) - G.*cos(theta));
                dscdS = zeros(1, length(S));
                
            otherwise
                error('MATLAB:contactodes:BadRegionIndex','Incorrect region index: %d',region);
        end
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS; dscdS];
        
    end

Odes = ContactEqns(x, M, region);

end
