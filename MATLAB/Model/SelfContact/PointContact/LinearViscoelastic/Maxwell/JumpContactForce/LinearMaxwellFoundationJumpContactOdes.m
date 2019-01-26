function Odes = LinearMaxwellFoundationJumpContactOdes(x, M, region, MOld, modelParams)

K = modelParams.K;
L = modelParams.L;
gamma = modelParams.gamma;
Es = modelParams.Es;
Eb = modelParams.Eb;
ext = modelParams.ext;
nu = modelParams.nu;
etaK = modelParams.etaK;
dt = modelParams.dt;

SOld = MOld.y(1,:);
XOld = MOld.y(2,:);
YOld = MOld.y(3,:);
mOld = MOld.y(7,:);
POld = modelParams.P;
uHatOld = modelParams.uHat;

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
        fC = M(9,:);

        % Interpolate the parameters if they are non-constant
        if (length(gamma) > 1)
            gammaInterp = interparc(x/2, MOld.x, gamma, 'linear');
            gamma = gammaInterp(:, 2)';
        end
        
        if (length(K) > 1)
            KInterp = interparc(x/2, MOld.x, K, 'linear');
            K = KInterp(:,2)';
        end
        
        if (length(Es) > 1)
            EsInterp = interparc(x/2, MOld.x, Es, 'linear');
            Es = EsInterp(:,2)';
        end
        
        if (length(Eb) > 1)
            EbInterp = interparc(x/2, MOld.x, Eb, 'linear');
            Eb = EbInterp(:,2)';
        end
        
        if (length(SOld) > 1)
            SOldInterp = interparc(x/2, MOld.x, SOld, 'linear');
            SOld = SOldInterp(:,2)';
        end
        
        if (length(XOld) > 1)
            XOldInterp = interparc(x/2, MOld.x, XOld, 'linear');
            XOld = XOldInterp(:,2)';        
        end
        
        if (length(YOld) > 1)
            YOldInterp = interparc(x/2, MOld.x, YOld, 'linear');
            YOld = YOldInterp(:,2)';        
        end
        
        if (length(mOld) > 1)
            mOldInterp = interparc(x/2, MOld.x, mOld, 'linear');
            mOld = mOldInterp(:,2)';        
        end
        
        if (length(POld) > 1)
            POldInterp = interparc(x/2, MOld.x, POld, 'linear');
            POld = POldInterp(:,2)';
        end
        
        if (length(uHatOld) > 1)
            uHatOldInterp = interparc(x/2, MOld.x, uHatOld, 'linear');
            uHatOld = uHatOldInterp(:,2)';
        end
        
        % If the model is extensible, set alpha to the tension, otherwise,
        % set alpha = 1.
        if (ext == 1)
            alpha = 1 + (F.*cos(theta) + G.*sin(theta))./Es;
        else
            alpha = 1;
        end
        
        dscdS = zeros(1, length(S));
        dfCdS = zeros(1, length(S));
        
        Delta = sqrt((X - S).^2 + (Y).^2);
        DeltaOld = sqrt((XOld - SOld).^2 + (YOld).^2);
        
        % Update the viscoelastic stress in time
        P = Delta - DeltaOld + (1 - dt*nu).*POld;
        
        % Update the intrinsic curvature
        uHat = uHatOld + etaK*dt.*mOld./Eb;
        
        switch region
            case 1 % ODEs for [0, sc1]
 
                dSdS = (sc).*ones(1, length(S));
                dxdS = (sc).*gamma.*alpha.*cos(theta);
                dydS = (sc).*gamma.*alpha.*sin(theta);
%                 dFdS = (sc).*K.*P.*(X - S);
%                 dGdS = (sc).*K.*P.*Y;

                dthetadS = (sc).*gamma.*(m./Eb + uHat);
                dmdS = (sc).*gamma.*alpha.*(F.*sin(theta) - G.*cos(theta));
                
            case 2 % ODEs for [sc2, L]
                
                dSdS = (L - sc).*ones(1, length(S));
                dxdS = (L - sc).*gamma.*alpha.*cos(theta);
                dydS = (L - sc).*gamma.*alpha.*sin(theta);
                dFdS = (L - sc).*K.*P.*(X - S);
                dGdS = (L - sc).*K.*P.*Y;              
                dthetadS = (L - sc).*gamma.*(m./Eb + uHat);
                dmdS = (L - sc).*gamma.*alpha.*(F.*sin(theta) - G.*cos(theta));
                
            otherwise
                error('MATLAB:contactodes:BadRegionIndex','Incorrect region index: %d',region);
        end
        
        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS; dscdS; dfCdS];
    end

Odes = ContactEqns(x, M, region);

end
