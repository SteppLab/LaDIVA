%%
% Rule2BodyCoverParameters: Function implementing the computation of the 
% geometrical and % biomechanical parameters of the Body Cover Model 
% resulting from the normalized activity levels for cricothyroid (CT), 
% thyroarytenoid (TA), and lateral cricoarytenoid (LC) muscles. This method
% allows for the dynamical simulation of muscle control of vocal folds 
% oscillations.
%
% Structure: BCMParam = Rule2BodyCoverParameters(BCMObj, a_CT, a_TA, a_LC), 
% where
%
% BCMObj: is an object from BodyCoverModel (handle) class, 
% a_CT (default 0.25): is the activation level for CT muscle, 
% a_TA (default 0.10): is the activation level for TA muscle,
% a_LC (default 0.50): is the activation level for LC muscle.
% BCMParam: Struct gathering the computed model parameters.
%
% References:
% [1] I. R. Titze and B. H. Story, “Rules for controlling low-dimensional 
%     vocal fold models with muscle activation,�? J. Acoust. Soc. Am., 
%     vol. 112, p. 1064, 2002.
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function BCMParam = Rule2BodyCoverParameters(BCMObj, varargin)
  % Default activation levels
  a_CT = 0.25; a_TA = 0.10; a_LC = 0.50;
  
  % Check the input: Normalized activation levels
  if (nargin>=2)
    if isnumeric(varargin{1})&&(varargin{1}>=0)&&(varargin{1}<=1)
      a_CT = varargin{1};
    else
      error('Activation level for CT muscle must be a real number in the range [0, 1]!')
    end
    
    if (nargin>=3)
      if isnumeric(varargin{2})&&(varargin{2}>=0)&&(varargin{2}<=1)
        a_TA = varargin{2};
      else
        error('Activation level for TA muscle must be a real number in the range [0, 1]!')
      end
      
      if (nargin==4)
        if isnumeric(varargin{3})&&(varargin{3}>=0)&&(varargin{3}<=1)
          a_LC = varargin{3};
        else
          error('Activation level for LC muscle must be a real number in the range [0, 1]!')
        end
      elseif  (nargin>4)
        error('More input variables than allowed. Read function description!')  
      end
    end
  end
      
  % Elongation rule [-]: Equation 55
  epsilon = MuscleActivation.G*(MuscleActivation.R*a_CT - a_TA) - MuscleActivation.H*a_LC;
  
  % Fold elongation : Equation 56
  Lg = BCMObj.Lg0*(1+epsilon); % (EN EL PAPER DICE L = L0*(1+epsilon) pero en la enmienda aparece L = L0/(1+epsilon);!!!!!) <------!!!!
%   Lg = BCMObj.Lg0/(1+epsilon); % (EN EL PAPER DICE L = L0*(1+epsilon) pero en la enmienda aparece L = L0/(1+epsilon);!!!!!) <------!!!! THIS IS WRONG!!
  
  % Thickness rule : Equation 59
  Tg = BCMObj.Tg0/(1+0.8*epsilon);
  
  % - Point rule : Equation 58
  Znodal = (1+a_TA)*Tg/3;
  
  % Depth rules
%   if strcmpi(BCMObj.sex, 'male')
%     Dmuc = MuscleActivation.DEPTH_MUC_MALE;
%     Dlig = MuscleActivation.DEPTH_LIG_MALE;
%     Dmus = MuscleActivation.DEPTH_MUS_MALE;
%   elseif strcmpi(BCMObj.sex, 'female')
%     Dmuc = MuscleActivation.DEPTH_MUC_FEMALE;
%     Dlig = MuscleActivation.DEPTH_LIG_FEMALE;
%     Dmus = MuscleActivation.DEPTH_MUS_FEMALE;
%   else
%     error('Error in sex especification in the Body Cover Model')  
%   end
  
  Dmuc = BCMObj.DMuc0;
  Dlig = BCMObj.DLig0;
  Dmus = BCMObj.DMus0;
  D_body = (a_TA*Dmus+0.5*Dlig)/(1+0.2*epsilon); % Equation 60
  D_cover = (Dmuc+0.5*Dlig)/(1+0.2*epsilon); % Euqation 61
  % Adduction rule : Equation 62
  xi_02 = 0.25*BCMObj.Lg0*(1-2.0*a_LC); % [m] Upper mass
  % Convergence rule [cm]
%     xi_c = T*(0.05-0.15*Ata); % Equation 63
%     % Nearly rectangular approach
  xi_c = Tg*tan(0.0001);
  xi_01 = xi_c + xi_02; % [m] Lower mass <- errata of Version 26.01.B (old = xi_c - xi_02);
  
  % Cover and Body Stress % Equation 53/ 54
  sigma_muc = PassiveStress(epsilon,MuscleActivation.EPSILON_1_MUC,MuscleActivation.EPSILON_2_MUC,...
                    MuscleActivation.SIGMA_0_MUC,MuscleActivation.SIGMA_2_MUC,MuscleActivation.C_MUC); % [Pa]
                
  sigma_lig = PassiveStress(epsilon,MuscleActivation.EPSILON_1_LIG,MuscleActivation.EPSILON_2_LIG,...
                    MuscleActivation.SIGMA_0_LIG,MuscleActivation.SIGMA_2_LIG,MuscleActivation.C_LIG); % [Pa]
                
  sigma_mus = a_TA*MuscleActivation.SIGMA_M_TA*max(0,1-MuscleActivation.B_TA*(epsilon-MuscleActivation.EPSILON_M_TA)^2) + ...
              PassiveStress(epsilon,MuscleActivation.EPSILON_1_TA,MuscleActivation.EPSILON_2_TA,...
                    MuscleActivation.SIGMA_0_TA,MuscleActivation.SIGMA_2_TA,MuscleActivation.C_TA); % [Pa]
  sigma_body = (0.5*sigma_lig*Dlig+sigma_mus*Dmus)/D_body; % [Pa] : Equation 51
  sigma_cover = (0.5*sigma_lig*Dlig+sigma_muc*Dmuc)/D_cover; % [Pa] ; Equation 52
  
  % Computation of Body cover model parameters 
  Tl = Znodal; % [m]
  Tu = Tg-Znodal; % [m]
  % Geometric and elastic parameters for the Cover (l: lower mass, u: upper mass, c: coupling)
  kl = 2*MuscleActivation.SHEARMODULUS_COVER*(Lg*Tg/D_cover)*(Znodal/Tg) ...
        + (pi^2)*sigma_cover*(D_cover/Lg)*Znodal; % [Pa*m] or [N/m] : Equation 44
  ku = 2*MuscleActivation.SHEARMODULUS_COVER*(Lg*Tg/D_cover)*(1-(Znodal/Tg)) ...
        + (pi^2)*sigma_cover*(D_cover/Lg)*Tg*(1-(Znodal/Tg)); % [Pa*m] or [N/m] : Equation 45
  kc = ((1/2)*MuscleActivation.SHEARMODULUS_COVER*(Lg*D_cover/Tg)*((1/3)-(Znodal/Tg)*(1-(Znodal/Tg)))^(-1) ...
        - 2*MuscleActivation.SHEARMODULUS_COVER*(Lg*Tg/D_cover))*(Znodal/Tg)*(1-(Znodal/Tg)); % [Pa*m] or [N/m] : Equation 46
  ml = MuscleActivation.TISSUE_DENS*Lg*Tg*D_cover*(Znodal/Tg); % [kg] : Equation 47
  mu = MuscleActivation.TISSUE_DENS*Lg*Tg*D_cover*(1-(Znodal/Tg)); % [kg] : Equation 48
  % Geometric and elastic parameters for the Body
  kb = 2*MuscleActivation.SHEARMODULUS_BODY*(Lg*Tg/D_body) ...
        + (pi^2)*sigma_body*(D_body/Lg)*Tg; % [Pa*m] or [N/m] : Equation 49
  mb = MuscleActivation.TISSUE_DENS*Lg*Tg*D_body; % [kg]: Equation 50
      
  % Updating the parameters of the Body Cover object
  BCMObj.a_CT = a_CT;
  BCMObj.a_TA = a_TA;
  BCMObj.a_LC = a_LC;
  BCMObj.epsilon = epsilon;
  BCMObj.Lg = Lg;
  BCMObj.Tg = Tg;
  BCMObj.Znodal = Znodal;
  BCMObj.Tl = Tl;
  BCMObj.Tu = Tu;
  BCMObj.kl = kl;
  BCMObj.ku = ku;
  BCMObj.kc = kc;
  BCMObj.ml = ml;
  BCMObj.mu = mu;
  BCMObj.kb = kb;
  BCMObj.mb = mb;
  BCMObj.xu0 = xi_02;
  BCMObj.xl0 = xi_01;
  BCMObj.xi_01 = xi_01;
  BCMObj.xi_02 = xi_02;
  BCMObj.ParamSet_MuscleRules = true;
  
  % Gathering the resulting parameters in the output struct
  BCMParam.epsilon = epsilon;
  BCMParam.Lg = Lg;
  BCMParam.Tg = Tg;
  BCMParam.Znodal = Znodal;
  BCMParam.Tl = Tl;
  BCMParam.Tu = Tu;
  BCMParam.kl = kl;
  BCMParam.ku = ku;
  BCMParam.kc = kc;
  BCMParam.ml = ml;
  BCMParam.mu = mu;
  BCMParam.kb = kb;
  BCMParam.mb = mb;
  BCMParam.xi_01 = xi_01;
  BCMParam.xi_02 = xi_02;
end

% Passive stress formula
function out = PassiveStress(epsilon,epsilon_1,epsilon_2,sigma_0,sigma_2,C)
  if epsilon < epsilon_1
    out = 0;
  elseif epsilon <= epsilon_2
    out = -(sigma_0/epsilon_1)*(epsilon-epsilon_1);
  else
    out = -(sigma_0/epsilon_1)*(epsilon-epsilon_1) + sigma_2*(exp(C*(epsilon-epsilon_2))-C*(epsilon-epsilon_2)-1);
  end
end