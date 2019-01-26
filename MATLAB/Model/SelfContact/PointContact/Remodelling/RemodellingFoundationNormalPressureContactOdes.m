function Odes = RemodellingFoundationNormalPressureContactOdes(x, M, region, MOld, modelParams)

K = modelParams.K;
L = modelParams.L;
gamma = modelParams.gamma;
Eb = modelParams.Eb;
nu = modelParams.nu;
dt = modelParams.dt;

XOld = MOld.y(2,:);
YOld = MOld.y(3,:);
PxOld = modelParams.Px;
PyOld = modelParams.Py;

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
        A = M(9,:);
        pC = M(10,:);

        % Interpolate the parameters if they are non-constant
        if (length(gamma) > 1)
            gammaInterp = interparc(x/4, MOld.x, gamma, 'linear');
            gamma = gammaInterp(:, 2)';
        end
        
        if (length(Eb) > 1)
            EbInterp = interparc(x/4, MOld.x, Eb, 'linear');
            Eb = EbInterp(:,2)';
        end

        if (length(XOld) > 1)
            XOldInterp = interparc(x/4, MOld.x, XOld, 'linear');
            XOld = XOldInterp(:,2)';        
        end
        
        if (length(YOld) > 1)
            YOldInterp = interparc(x/4, MOld.x, YOld, 'linear');
            YOld = YOldInterp(:,2)';        
        end
        
        if (length(PxOld) > 1)
            PxOldInterp = interparc(x/4, MOld.x, PxOld, 'linear');
            PxOld = PxOldInterp(:,2)';
        end
        
        if (length(PyOld) > 1)
            PyOldInterp = interparc(x/4, MOld.x, PyOld, 'linear');
            PyOld = PyOldInterp(:,2)';
        end
        
        dscdS = zeros(1, length(S));
        dpCdS = zeros(1, length(S));
        
        % Update the foundation shape in time
        Px = PxOld + dt*nu.*(XOld - PxOld);
        Py = PyOld + dt*nu.*(YOld - PyOld); 
       
        
        switch region
            case 1 % ODEs for [0, sc1]
                
                l = sc;
 
                dSdS = (l).*ones(1, length(S));
                dxdS = (l).*gamma.*cos(theta);
                dydS = (l).*gamma.*sin(theta);
                dFdS = (l).*gamma.*K.*(X - Px);
                dGdS = (l).*gamma.*K.*(Y - Py);
                dthetadS = (l).*gamma.*(m./Eb);
                dmdS = (l).*gamma.*(F.*sin(theta) - G.*cos(theta));
                dAdS = zeros(1, length(S));

            case 2 % ODEs for [sc1, L/2]
                
                l = 0.5*L - sc;
                
                dSdS = (l).*ones(1, length(S));
                dxdS = (l).*gamma.*cos(theta);
                dydS = (l).*gamma.*sin(theta);
                dFdS = (l).*gamma.*(K.*(X - Px) + pC.*sin(theta));
                dGdS = (l).*gamma.*(K.*(Y - Py) - pC.*cos(theta));                           
                dthetadS = (l).*gamma.*(m./Eb);
                dmdS = (l).*gamma.*(F.*sin(theta) - G.*cos(theta));
                dAdS = (l).*gamma.*X.*sin(theta);

            case 3 % ODEs for [L/2, sc2]
                
                l = 0.5*L - sc;
                
                dSdS = (l).*ones(1, length(S));
                dxdS = (l).*gamma.*cos(theta);
                dydS = (l).*gamma.*sin(theta);
                dFdS = (l).*gamma.*(K.*(X - Px) + pC.*sin(theta));
                dGdS = (l).*gamma.*(K.*(Y - Py) - pC.*cos(theta));                           
                dthetadS = (l).*gamma.*(m./Eb);
                dmdS = (l).*gamma.*(F.*sin(theta) - G.*cos(theta));
                dAdS = (l).*gamma.*X.*sin(theta);
                
            case 4 % ODEs for [sc2, L]
                
                l = sc;
                
                dSdS = (l).*ones(1, length(S));
                dxdS = (l).*gamma.*cos(theta);
                dydS = (l).*gamma.*sin(theta);
                dFdS = (l).*gamma.*K.*(X - Px);
                dGdS = (l).*gamma.*K.*(Y - Py);                           
                dthetadS = (l).*gamma.*(m./Eb);
                dmdS = (l).*gamma.*(F.*sin(theta) - G.*cos(theta));
                dAdS = zeros(1, length(S));
                
            otherwise
                error('MATLAB:contactodes:BadRegionIndex','Incorrect region index: %d',region);
        end
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS; dscdS; dAdS; dpCdS];
    end

Odes = ContactEqns(x, M, region);

end
