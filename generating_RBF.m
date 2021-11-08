%% Generating adial basis function for lookuptable
%  input lookuptable load('BCMlookuptable.mat')
%  newrbe(X,T)
%  X = [RxQ] - Q input vectors
%  T = [SxQ] - S r=target class vectors
% SPREAD = spread of rbf - default is 1.0 (no errors)

  % --- Estimating the model ---
%         % INPUT:
%         % Motor_new : [Nd,Ns] list of Ns samples of motor configurations
%         % Motor_LF  : [3,Ns] list of Ns samples of LF model parameters associated with the same Motor_new samples
%         ncl=32;               % number of centroids
%         maxiter=1000;         % maximum number of iterations
%         gmfit_art=gmdistribution.fit(Motor_new',ncl,'Start','randSample','Replicates',1,'SharedCov',true,'Regularize',1e-5,'Options',struct('Display','iter','MaxIter',maxiter)); % estimates model
%         MODEL=struct('mu',gmfit_art.mu,'iSigma',pinv(gmfit_art.Sigma),'p',gmfit_art.PComponents.','beta',gmfit_art.posterior(Motor_new')\Motor_LF');
%
%         % --- Using an estimated model to approximate the original transformation ---
%         % INPUT: motor_new : Nd element vector of motor configurations
%         % OUTPUT: motor_LF  : 3 element vector of fitted LF model parameters
%         dx=bsxfun(@minus,motor_new(:)',MODEL.mu);
%         p=-sum((dx*MODEL.iSigma).*dx,2)/2;
%         p=MODEL.p.*exp(p-max(p));
%         p=p/sum(p);
%         motor_LF=MODEL.beta'*p;

Motor_new(1,:) = lookuptable.aCTrange;  % 0:0.001:1 %1001 points
Motor_new(2,:) = lookuptable.aTArange;  % 0:0.001:1 %1001 points
Motor_new(3,:) = lookuptable.Psrange;   % 10:2:2010 %1001 points

Motor_LF(1,:) = lookuptable.F0; %Hz
Motor_LF(2,:) = lookuptable.SPL; %dB
Motor_LF(3,:) = ones(1,lenth(Motor_LF)); %voicing - always 1


%net_F0 = newrbe(X,T_F0);
%net_SPL = newrbe(X,T_SPL);

