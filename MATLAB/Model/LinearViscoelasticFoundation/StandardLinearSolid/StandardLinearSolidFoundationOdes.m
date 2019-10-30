function Odes = StandardLinearSolidFoundationOdes(x, M, MOld, parameters)

y0 = parameters.y0;
K = parameters.K;
beta = parameters.beta;
L = parameters.L;
nu = parameters.nu;
dt = parameters.dt;
gamma = parameters.gamma;
Eb = parameters.Eb;

SOld = MOld.y(1,:);
XOld = MOld.y(2,:);
YOld = MOld.y(3,:);
POld = parameters.P;

    function dMdS = VicoelasticFoundationEqns(x, M)
        
        % Define the state variables
        S = M(1,:);
        X = M(2,:);
        Y = M(3,:);
        F = M(4,:);
        G = M(5,:);
        theta = M(6,:);
        m = M(7,:);
        
        % Interpolate the parameters if they are non-constant
        if (length(gamma) > 1)
            gamma = interp1(MOld.x, gamma, x);
        end
        
        if (length(Eb) > 1)
            Eb = interp1(MOld.x, Eb, x);
        end
        
        if (length(SOld) > 1)
            SOld = interp1(MOld.x, SOld, x);
        end
        
        if (length(XOld) > 1)
            XOld = interp1(MOld.x, XOld, x);
        end
        
        if (length(YOld) > 1)
            YOld = interp1(MOld.x, YOld, x);
        end
        
        if (length(POld) > 1)
            POld = interp1(MOld.x, POld, x);
        end
        
        DeltaOld = sqrt((XOld - SOld).^2 + (YOld).^2);
        Delta = sqrt((X - S).^2 + (Y).^2);
        
        P = (1 + dt*beta*nu).*(POld) + dt*beta*(1 - beta)*K*nu.*(DeltaOld - y0) ...
        + (K).*(Delta - DeltaOld);
        
        dSdS = L.*ones(1, length(S));
        dxdS = L.*gamma.*cos(theta);
        dydS = L.*gamma.*sin(theta);
        dFdS = L.*gamma.*P.*(X - S)./Delta;
        dGdS = L.*gamma.*P.*(Y)./Delta;
        dthetadS = L.*gamma.*(m./Eb);
        dmdS = L.*gamma.*(F.*sin(theta) - G.*cos(theta));
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
        
    end

Odes = VicoelasticFoundationEqns(x, M);

end
