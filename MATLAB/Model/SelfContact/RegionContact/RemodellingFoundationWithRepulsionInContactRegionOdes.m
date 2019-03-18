function Odes = RemodellingFoundationWithRepulsionInContactRegionOdes(x, M, MOld, modelParams)

K = modelParams.K;
gamma = modelParams.gamma;
Eb = modelParams.Eb;
Q = modelParams.Q;
N = modelParams.N;
L = modelParams.L;

xOld = MOld.x;

Px = modelParams.Px;
Py = modelParams.Py;

    function dMdS = ContactEqns(x, M)
        
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
        
        
        l = sc;
        
        dSdS = (l).*ones(1, length(S));
        dxdS = (l).*gamma.*cos(theta);
        dydS = (l).*gamma.*sin(theta);
        dFdS = (l).*gamma.*(K.*(X - Px) + Q./(X - 1.01*L).^N.*(X < L));
        dGdS = (l).*gamma.*K.*(Y - Py);
        dthetadS = (l).*gamma.*(m./Eb);
        dmdS = (l).*gamma.*(F.*sin(theta) - G.*cos(theta));
        dscdS = zeros(1, length(S));
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS; dscdS];
        
    end

Odes = ContactEqns(x, M);

end
