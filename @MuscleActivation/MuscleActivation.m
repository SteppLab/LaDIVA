classdef MuscleActivation
% Class implementing the muscle rules for computing the mechanical
% parameters controlling the oscillations of body cover model of the vocal
% folds. Current class version include static methods for computing the 
% geometric and mechanical parameters involved in body cover formulation
% as a function of the normalized activatity levels for 
% cricothyroid (CT), thyroarytenoid (TA), and lateral cricoarytenoid (LC)
% muscles.
%
% References:
% [1] I. R. Titze and B. H. Story, “Rules for controlling low-dimensional 
%     vocal fold models with muscle activation,” J. Acoust. Soc. Am., 
%     vol. 112, p. 1064, 2002.
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
%
  properties (Constant, Hidden)
    % Unisex Constants
    G = 0.2; % Gain of elongation [-]
    R = 3.0; % Torque ratio [-]
    H = 0.2; % Adductory strain factor [-]
    TISSUE_DENS = 1040; % [kg/m^3] Tissue density
    SHEARMODULUS_COVER = 500; % [Pa] Shear modulus of cover layers (0.500 [kPa])
    SHEARMODULUS_BODY = 100; % [Pa] Shear modulus of body layer (1 [kPa])
%     % Sex dependent Constants
%     DEPTH_MUC_MALE = 0.2e-2; % [m] Depth of mucosa (0.2 [cm] in males and 0.15 [cm] in females)
%     DEPTH_LIG_MALE = 0.2e-2; % [m] Depth of ligament (0.2 [cm] in males and 0.15 [cm] in females)
%     DEPTH_MUS_MALE = 0.4e-2; % [m] Depth of TA muscle (0.4 [cm] in males and 0.3 [cm] in females)
%     DEPTH_MUC_FEMALE = 0.15e-2; % [m] Depth of mucosa (0.2 [cm] in males and 0.15 [cm] in females)
%     DEPTH_LIG_FEMALE = 0.15e-2; % [m] Depth of ligament (0.2 [cm] in males and 0.15 [cm] in females)
%     DEPTH_MUS_FEMALE = 0.3e-2; % [m] Depth of TA muscle (0.4 [cm] in males and 0.3 [cm] in females)
    % Stress-strain Constants
    % Mucosa
    EPSILON_1_MUC = -0.5; % [-]
    EPSILON_2_MUC = -0.35; % [-] (EN EL PAPER DICE +0,35!!!!!) <------!!!!
    SIGMA_0_MUC   = 500; % [Pa]       0.5 [kPa]
    SIGMA_2_MUC   = 30000; % [Pa]     30 [kPa]
    C_MUC         = 4.4; % [-]
    % Ligament
    EPSILON_1_LIG = -0.5; % [-]
    EPSILON_2_LIG = -0.00; % [-]
    SIGMA_0_LIG   = 400; % [Pa]      0.4 [kPa]
    SIGMA_2_LIG   = 1393; % [Pa]     1.393 [kPa] 
    C_LIG         = 17.0; % [-]
    % TA Muscle
    EPSILON_1_TA = -0.5; % [-]
    EPSILON_2_TA = -0.05; % [-]
    SIGMA_0_TA   = 1000; % [Pa]     1.0 [kPa]
    SIGMA_2_TA   = 1500; % [Pa]     1.5 [kPa]
    C_TA         = 6.5; % [-]
    SIGMA_M_TA   = 105000; % [Pa] maximum active stress in the TA muscle 105 [kPa]
    EPSILON_M_TA = 0.4; % [-]
    B_TA         = 1.07; % [-]
  end
  methods (Static)
    function Param = getBodyCoverModelParameters(BCMObj)
    % Function for showing the Body Cover parameters handled by the
    % MuscleActivation methods.
      Param=struct;
      try
        Param.Lg = BCMObj.Lg; % [m] Dynamic vocal fold length
        Param.Tg = BCMObj.Tg; % [m] Dynamic vocal fold thickness
        Param.Tu = BCMObj.Tu; % [m] Dynamic thickness for the upper mass
        Param.Tl = BCMObj.Tl; % [m] Dynamic thickness for the lower mass
        Param.epsilon = BCMObj.epsilon; % [-] Vocal fold elongation
        Param.Znodal = BCMObj.Znodal;  % [m] Nodal point on the medial surface
        Param.mu = BCMObj.mu; % [kg] Mass of the upper cover mass
        Param.ml = BCMObj.ml; % [kg] Mass of the lower cover mass
        Param.mb = BCMObj.mb; % [kg] Mass of the body mass
        Param.ku = BCMObj.ku; % [N/m] Linear spring constant for the upper mass
        Param.kl = BCMObj.kl; % [N/m] Linear spring constant for the lower mass
        Param.kc = BCMObj.kc; % [N/m] Coupling spring constant for the cover layer
        Param.kb = BCMObj.kb; % [N/m] Linear spring constant for the body mass
        Param.a_CT = BCMObj.a_CT; % [-] Activity level for cricothyroid (CT) muscle
        Param.a_TA = BCMObj.a_TA; % [-] Activity level for thyroarytenoid (TA) muscle
        Param.a_LC = BCMObj.a_LC; % [-] Activity level for lateral cricoarytenoid (LC) muscle
      catch
        error('The input is not a BodyCoverModel object!')
      end
    end
    
    BCMParam = Rule2BodyCoverParameters(BCMObj, varargin)
    
  end    
end

