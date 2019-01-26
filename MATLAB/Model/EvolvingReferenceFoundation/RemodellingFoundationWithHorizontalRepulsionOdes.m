function Odes = RemodellingFoundationWithHorizontalRepulsionOdes(x, M, MOld, parameters)

K = parameters.K;
L = parameters.L;
Q = parameters.Q;
N = parameters.N;
gamma = parameters.gamma;
Eb = parameters.Eb;

Px = parameters.Px;
Py = parameters.Py;

    function dMdS = FoundationEqns(x, M)
        
        % Define the state variables
        S = M(1,:);
        X = M(2,:);
        Y = M(3,:);
        F = M(4,:);
        G = M(5,:);
        theta = M(6,:);
        m = M(7,:);
        
        if (length(gamma) > 1)
            gamma = interp1(MOld.x, gamma, x);
        end
        
        if (length(K) > 1)
            K = interp1(MOld.x, K, x);
        end
        
        if (length(Eb) > 1)
            Eb = interp1(MOld.x, Eb, x);
        end
        
        if (length(Px) > 1)
            Px = interp1(MOld.x, Px, x);
        end
        
        if (length(Py) > 1)
            Py = interp1(MOld.x, Py, x);
        end
        
        dSdS = L.*ones(1, length(S));
        dxdS = L.*gamma.*cos(theta);
        dydS = L.*gamma.*sin(theta);
        dFdS = L.*gamma.*(K.*(X - Px) + Q./(X - 1.01*L).^N.*(X < L));
        dGdS = L.*gamma.*K.*(Y - Py);
        dthetadS = L.*gamma.*(m./Eb);
        dmdS = L.*gamma.*(F.*sin(theta) - G.*cos(theta));

        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
        
    end

Odes = FoundationEqns(x, M);

end
