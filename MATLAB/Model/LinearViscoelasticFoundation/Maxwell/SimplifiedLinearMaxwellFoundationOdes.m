function Odes = SimplifiedLinearMaxwellFoundationOdes(x, M, MOld, parameters)

K = parameters.K;
L = parameters.L;
nu = parameters.nu;
etaK = parameters.etaK;
dt = parameters.dt;
gamma = parameters.gamma;
Es = parameters.Es;
Eb = parameters.Eb;
ext = parameters.ext;

YOld = MOld.y(3,:);
mOld = MOld.y(7,:);
POld = parameters.P;
uHatOld = parameters.uHat;

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
            gamma = interp1(MOld.x, gamma, x).*(length(gamma) > 1);
        end
        
        if (length(gamma) > 1)
            gamma = interp1(MOld.x, gamma, x);
        end
        
        if (length(K) > 1)
            K = interp1(MOld.x, K, x);
        end
        
        if (length(Es) > 1)
            Es = interp1(MOld.x, Es, x);
        end
        
        if (length(Eb) > 1)
            Eb = interp1(MOld.x, Eb, x);
        end
        
        if (length(YOld) > 1)
            YOld = interp1(MOld.x, YOld, x);
        end
        
        if (length(mOld) > 1)
            mOld = interp1(MOld.x, mOld, x);
        end
        
        if (length(POld) > 1)
            POld = interp1(MOld.x, POld, x);
        end
       
        if (length(uHatOld) > 1)
            uHatOld = interp1(MOld.x, uHatOld, x);
        end
       
        
        % If the model is extensible, set alpha to the tension, otherwise,
        % set alpha = 1.
        if (ext == 1)
            alpha = 1 + (F.*cos(theta) + G.*sin(theta))./Es;
        else
            alpha = 1;
        end
        
        % Update the viscoelastic stress in time
        P = Y - YOld + (1 - dt*nu).*POld;
        
        % Update the intrinsic curvature
        uHat = uHatOld + etaK*dt.*mOld./Eb;
        
        dSdS = L.*ones(1, length(S));
        dxdS = L.*gamma.*alpha.*cos(theta);
        dydS = L.*gamma.*alpha.*sin(theta);
        dFdS = 0;
        dGdS = L*K.*P;
        dthetadS = L.*gamma.*(m./Eb + uHat);
        dmdS = L.*gamma.*alpha.*(F.*sin(theta) - G.*cos(theta));

        dMdS = [dSdS; dxdS; dydS; dFdS; dGdS; dthetadS; dmdS];
        
    end

Odes = VicoelasticFoundationEqns(x, M);

end
