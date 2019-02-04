function Odes = RemodellingFoundationNormalPressureContactRegionOdes(x, M, region, MOld, modelParams)

K = modelParams.K;
gamma = modelParams.gamma;
Eb = modelParams.Eb;

contIndex = find(MOld.x == 1);
xOld = MOld.x([1:contIndex, (contIndex + 2):end]);

Px = modelParams.Px;
Py = modelParams.Py;

    function dMdS = ContactEqns(x, M, region)
        
        % Define the state variables
%         S = M(1,:);
        X = M(1,:);
        Y = M(2,:);
        F = M(3,:);
        G = M(4,:);
        theta = M(5,:);
        m = M(6,:);
        sc = M(7,:);
        A = M(8,:);
        pC = M(9,:);
        
        % Interpolate the parameters if they are non-constant
        if (length(gamma) > 1)
            gammaOld = gamma([1:contIndex, (contIndex + 2):end]);
            gamma = interp1(xOld, gammaOld, x);
        end
        
        if (length(Eb) > 1)
            EbOld = Eb([1:contIndex, (contIndex + 2):end]);
            Eb = interp1(xOld, EbOld, x);
        end
        
        if (length(Px) > 1)
            PxOld = Px([1:contIndex, (contIndex + 2):end]);
            Px = interp1(xOld, PxOld, x);
        end
        
        if (length(Py) > 1)
            PyOld = Py([1:contIndex, (contIndex + 2):end]);
            Py = interp1(xOld, PyOld, x);
        end

        switch region
            case 1 % ODEs for [0, sc1]
                
                l = sc;
 
                dxdS = (l).*gamma.*cos(theta);
                dydS = (l).*gamma.*sin(theta);
                dFdS = (l).*gamma.*(K.*(X - Px) - pC.*sin(theta));
                dGdS = (l).*gamma.*(K.*(Y - Py) + pC.*cos(theta));
                dthetadS = (l).*gamma.*(m./Eb);
                dmdS = (l).*gamma.*(F.*sin(theta) - G.*cos(theta));
                dscdS = zeros(1, length(S));
                dAdS = -(l).*X.*gamma.*sin(theta);
                dpCdS = zeros(1, length(S));

            case 2 % ODEs for [sc1, L/2]
                
                l = 0.5 - sc;
                
                dxdS = (l).*gamma.*cos(theta);
                dydS = (l).*gamma.*sin(theta);
                dFdS = (l).*gamma.*K.*(X - Px);
                dGdS = (l).*gamma.*K.*(Y - Py);                           
                dthetadS = (l).*gamma.*(m./Eb);
                dmdS = (l).*gamma.*(F.*sin(theta) - G.*cos(theta));
                dscdS = zeros(1, length(S));
                dAdS = zeros(1, length(S));
                dpCdS = zeros(1, length(S));

            otherwise
                error('MATLAB:contactodes:BadRegionIndex','Incorrect region index: %d',region);
        end
        
        dMdS = [dxdS; dydS; dFdS; dGdS; dthetadS; dmdS; dscdS; dAdS; dpCdS];
        
    end

Odes = ContactEqns(x, M, region);

end
