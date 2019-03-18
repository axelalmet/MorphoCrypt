function Odes = RemodellingFoundationNormalPressureOutContactRegionOdes(x, M, MOld, modelParams)

K = modelParams.K;
gamma = modelParams.gamma;
Eb = modelParams.Eb;

xOld = MOld.x;

Px = modelParams.Px;
Py = modelParams.Py;

sc = modelParams.sc;

    function dMdS = ContactEqns(x, M)
        
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
            gamma = interp1(xOld, gamma, x);
        end
        
        if (length(Eb) > 1)
            Eb = interp1(xOld, Eb, x);
        end
        
        if (length(Px) > 1)
            Px = interp1(xOld, Px, x);
        end
        
        if (length(Py) > 1)
            Py = interp1(xOld, Py, x);
        end
        
        l = 0.5 - sc;
        
        dSdS = (l).*ones(1, length(S));
        dxdS = (l).*gamma.*cos(theta);
        dydS = (l).*gamma.*sin(theta);
        dFdS = (l).*gamma.*(K.*(X - Px));
        dGdS = (l).*gamma.*(K.*(Y - Py));
        dthetadS = (l).*gamma.*(m./Eb);
        dmdS = (l).*gamma.*(F.*sin(theta) - G.*cos(theta));
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
        
    end

Odes = ContactEqns(x, M);

end
