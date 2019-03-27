function Odes = RemodellingFoundationWithRepulsionSelfContactRegionOdes(x, M, region, MOld, modelParams)

K = modelParams.K;
gamma = modelParams.gamma;
Eb = modelParams.Eb;
Q = modelParams.Q;
N = modelParams.N;
L = modelParams.L;

splitIndex = find(MOld.x == 1, 1);
xOld = MOld.x([1:splitIndex, (splitIndex + 2):end]);

Px = modelParams.Px;
Py = modelParams.Py;

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
        sr = M(9,:);
        
        % Interpolate the parameters if they are non-constant
        if (length(gamma) > 1)
            gammaOld = gamma([1:splitIndex, (splitIndex + 2):end]);
            gamma = interp1(xOld, gammaOld, x);
        end
        
        if (length(Eb) > 1)
            EbOld = Eb([1:splitIndex, (splitIndex + 2):end]);
            Eb = interp1(xOld, EbOld, x);
        end
        
        if (length(Px) > 1)
            PxOld = Px([1:splitIndex, (splitIndex + 2):end]);
            Px = interp1(xOld, PxOld, x);
        end
        
        if (length(Py) > 1)
            PyOld = Py([1:splitIndex, (splitIndex + 2):end]);
            Py = interp1(xOld, PyOld, x);
        end
        
        dscdS = zeros(1, length(S));
        dsrdS = zeros(1, length(S));
        
        switch region
            case 1
                
                l = sc;
                
                dSdS = l.*ones(1, length(S));
                dxdS = l.*gamma.*cos(theta);
                dydS = l.*gamma.*sin(theta);
                dFdS = l.*gamma.*(K.*(X - Px) + Q./(X - 1.01*L).^N.*(X < L));
                dGdS = l.*gamma.*K.*(Y - Py);
                dthetadS = l.*gamma.*(m./Eb);
                dmdS = l.*gamma.*(F.*sin(theta) - G.*cos(theta));
                
            case 2
                
                l = L - (sc + sr);
                
                dSdS = l.*ones(1, length(S));
                dxdS = l.*gamma.*cos(theta);
                dydS = l.*gamma.*sin(theta);
                dFdS = l.*gamma.*(K.*(X - Px) + Q./(X - 1.01*L).^N.*(X < L));
                dGdS = l.*gamma.*K.*(Y - Py);
                dthetadS = l.*gamma.*(m./Eb);
                dmdS = l.*gamma.*(F.*sin(theta) - G.*cos(theta));
                
                
            otherwise
                error('MATLAB:contactodes:BadRegionIndex','Incorrect region index: %d',region);
        end
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS; dscdS; dsrdS];
        
    end

Odes = ContactEqns(x, M, region);

end