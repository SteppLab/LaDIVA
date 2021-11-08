classdef BodyCoverModel < handle
% Handle class for modeling vocal fold oscillations according to the
% lump-elements body-cover model of the vocal folds. The model describes
% the cover layer through two (lower and upper) masses, and the body part
% with a deep mass. All the three masses are conected each other. 
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body�?cover 
%     model of the vocal folds,�? J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
%
% Coded by Gabriel Alzamendi, July 2019 (Modified December 2019).
% Based on previous code by Matías Zañartu and Gabriel Galindo.
  
  properties (Constant, Hidden)
    % Constant parameters describign the anatomical average dimension for
    % normal male and female vocal folds.
    LREST_MALE = 1.6e-2; % [m] Vocal fold length at rest for male subject
    TREST_MALE = 0.3e-2; % [m] Vocal fold thickness at rest for male subject
    DEPTH_MUC_MALE = 0.2e-2; % [m] Depth of mucosa (0.2 [cm] in males and 0.15 [cm] in females)
    DEPTH_LIG_MALE = 0.2e-2; % [m] Depth of ligament (0.2 [cm] in males and 0.15 [cm] in females)
    DEPTH_MUS_MALE = 0.4e-2; % [m] Depth of TA muscle (0.4 [cm] in males and 0.3 [cm] in females)
      
    LREST_FEMALE = 1.0e-2; % [m] Vocal fold length at rest for female subject
    TREST_FEMALE = 0.2e-2; % [m] Vocal fold thickness at rest for female subject
    DEPTH_MUC_FEMALE = 0.15e-2; % [m] Depth of mucosa (0.2 [cm] in males and 0.15 [cm] in females)
    DEPTH_LIG_FEMALE = 0.15e-2; % [m] Depth of ligament (0.2 [cm] in males and 0.15 [cm] in females)
    DEPTH_MUS_FEMALE = 0.3e-2; % [m] Depth of TA muscle (0.4 [cm] in males and 0.3 [cm] in females)
  end
  
  properties (SetAccess = protected)
    % Dynamic state variable
    Model = 'BCM';
    xData = zeros(9,1); % State variable descring vocal fold oscillations:
                   % xData = [xu xl xb vu vl vb au al ab], where
                   %   x_: mass displacement in [m]
                   %   v_: mass velocity in [m/s]
                   %   a_: mass acceleration in [m/s^2]
    n_IterCont = 0; % Simulation time index
    % General object description
    Lg0 = 0; % [m] Vocal fold length at rest
    Tg0 = 0;% [m] Vocal fold thickness at rest
    sex = '';
    scale = 1.0; % Scale factor for converting the volume of the vocal fold (VFvol_out = scale * VFvol_rest)
    
    % Simulation parameters
    fs = 0; % [Hz] Sampling frequency for the numerical simulation
    Ts = 0; % [s] Sampling period for the numerical simulation
    SimParamOK = false;  % If 'true' simulation parameters are set,
                         % otherwise simulation parameters are missing
    NonLinMode = 1; % Binari factor setting on (=1) or off (=0) the 
                    % non-linear terms in the mechanical rules.
  end
    
  properties (GetAccess = protected)
    % Structural parameters
    Lg_init = 0; % [m] Vocal fold length at rest
    Tg_init = 0;% [m] Vocal fold thickness at rest
    DMuc_init = 0;% [m] Mucosa depth at rest
    DLig_init = 0;% [m] Ligament depth at rest
    DMus_init = 0;% [m] TA muscle depth at rest
    
    % Biomechanical parameters
    etau = 100*1e4; % [1/m^2] Non-linear spring coefficient for the upper mass
    etal = 100*1e4; % [1/m^2] Non-linear spring coefficient for the lower mass
    etab = 100*1e4; % [1/m^2] Non-linear spring coefficient for the body mass
    zetau0 = 0.6; % 0.4; % [-] Basic damping factor for the upper mass
    zetal0 = 0.1; % 0.4; % [-] Basic damping factor for the lower mass
    zetab = 0.2; % 0.2; % [-] Damping factor for the body mass
    
    % Collision biomechanical parameters
    xu_col = 0; % [m] Displacement where collision occurs for upper mass
    xl_col = 0; % [m] Displacement where collision occurs for lower mass
    etau_col = 500*1e4; % [1/m^2] Non-linear spring coefficient during collision for the upper mass
    etal_col = 500*1e4; % [1/m^2] Non-linear spring coefficient during collision for the lower mass
    
    % Boolean variables controlling BodyCoverModel object behavior
    symmetric = true; % If 'true' symmetric lateral displacements are considered, 
                      %  otherwise asymmetric displacement are assumed (INCOMPLETE)
    UseUpdated_ContactRules = false; % If 'true' new contact rules are used,
                                     %  otherwise original expressions are applied
    UseUpdated_AeroDrivingForces = false; % If 'true' new rules for computing the aerodynamic driving pressure are used, 
                                          % otherwise original expression are applied
  end
  
  properties (GetAccess = public, SetAccess = ?MuscleActivation)
    % Structural parameters  
    Lg = 0; % [m] Dynamic vocal fold length
    Tg = 0; % [m] Dynamic vocal fold thickness
    Tu = 0; % [m] Dynamic thickness for the upper mass
    Tl = 0; % [m] Dynamic thickness for the lower mass
    DMuc0 = 0;% [m] Mucosa depth at rest
    DLig0 = 0;% [m] Ligament depth at rest
    DMus0 = 0;% [m] TA muscle depth at rest
    epsilon = 0; % [-] Vocal fold elongation
    Znodal = 0; % [m] Nodal point on the medial surface
        
    % Normalized activity levels for the laryngeal muscles controlling the 
    % oscillations of the body Cover model
    a_CT = 0; % [-] Activity level for cricothyroid (CT) muscle
    a_TA = 0; % [-] Activity level for thyroarytenoid (TA) muscle
    a_LC = 0.5; % [-] Activity level for lateral cricoarytenoid (LC) muscle
    
    % Biomechanical parameters : Table 1 
    mu = 0; % [kg] Mass of the upper cover mass
    ml = 0; % [kg] Mass of the lower cover mass
    mb = 0; % [kg] Mass of the body mass
    xu0 = 0; % [m] Initial position of upper mass
    xl0 = 0; % [m] Initial position of lower mass
    xb0 = 3e-3; % [m] Initial position of body mass
    ku = 0; % [N/m] Linear spring constant for the upper mass
    kl = 0; % [N/m] Linear spring constant for the lower mass
    kc = 0; % [N/m] Coupling spring constant for the cover layer
    kb = 0; % [N/m] Linear spring constant for the body mass
    xi_01 = 0; % [m] Vocal process horizontal displacement for the lower mass 
    xi_02 = 0; % [m] Vocal process horizontal displacement for the upper mass
    % Boolean variables controlling BodyCoverModel object behavior
    ParamSet_MuscleRules = false;  % If 'true' model parameters are set trough muscle rules, 
                                   % otherwise model parameters need to be computed
  end
  
  properties (Dependent)%, Access = private)
    % Biomechanical parameters
    du % [Ns/m] Damping coefficient for the upper mass
    dl % [Ns/m] Damping coefficient for the lower mass
    db % [Ns/m] Damping coefficient for the body mass
    zetau % [-] Damping factor for the upper mass
    zetal % [-] Damping factor for the lower mass
    
    % Collision biomechanical parameters
    hu_col % [N/m] Linear spring constant during collision for the upper mass
    hl_col % [N/m] Linear spring constant during collision for the upper mass
    
    % Area variables
    au % [m2] Glottal area for the upper portion
    al % [m2] Glottal area for the lower portion
    ag % [m2] Proyected glottal area
    acont % [m2] Glottal contact area
  end
  
  methods
    % Class constructor
    function BCMObj = BodyCoverModel(varargin)
      if nargin ==0
        BCMObj.sex = 'male';
      elseif (nargin == 1)&&(ischar(varargin{1}))
        if (strcmpi(varargin{1},'male'))
          BCMObj.sex = 'male';
          BCMObj.Lg_init = BodyCoverModel.LREST_MALE;
          BCMObj.Tg_init = BodyCoverModel.TREST_MALE;
          BCMObj.DMuc_init = BodyCoverModel.DEPTH_MUC_MALE;
          BCMObj.DLig_init = BodyCoverModel.DEPTH_LIG_MALE;
          BCMObj.DMus_init = BodyCoverModel.DEPTH_MUS_MALE;
          InitModel(BCMObj);
          scalingVocalFold(BCMObj);
        elseif (strcmpi(varargin{1},'female'))
          BCMObj.sex = 'female';
          BCMObj.Lg_init = BodyCoverModel.LREST_FEMALE;
          BCMObj.Tg_init = BodyCoverModel.TREST_FEMALE;
          BCMObj.DMuc_init = BodyCoverModel.DEPTH_MUC_FEMALE;
          BCMObj.DLig_init = BodyCoverModel.DEPTH_LIG_FEMALE;
          BCMObj.DMus_init = BodyCoverModel.DEPTH_MUS_FEMALE;
          InitModel(BCMObj);
          scalingVocalFold(BCMObj);
        else
          error('Acceptable ''Sex'' especification are ''male'' or ''female'' ')
        end
      elseif (nargin == 1)&&(isa(varargin{1},'BodyCoverModel'))
        BCMAux=BodyCoverModel(varargin{1}.sex);
        BCMObj.sex = varargin{1}.sex;
        BCMObj.Lg_init = BCMAux.Lg_init;
        BCMObj.Tg_init = BCMAux.Tg_init;
        BCMObj.DMuc_init = BCMAux.DMuc_init;
        BCMObj.DLig_init = BCMAux.DLig_init;
        BCMObj.DMus_init = BCMAux.DMus_init;
        BCMObj.xData = varargin{1}.xData;
        BCMObj.n_IterCont = varargin{1}.n_IterCont;
        BCMObj.Lg0 = varargin{1}.Lg0;
        BCMObj.Tg0 = varargin{1}.Tg0;
        BCMObj.scale = varargin{1}.scale;
        BCMObj.fs = varargin{1}.fs;
        BCMObj.Ts = varargin{1}.Ts;
        BCMObj.SimParamOK = varargin{1}.SimParamOK;
        BCMObj.Lg = varargin{1}.Lg;
        BCMObj.Tg = varargin{1}.Tg;
        BCMObj.Tu = varargin{1}.Tu;
        BCMObj.Tl = varargin{1}.Tl;
        BCMObj.DMuc0 = varargin{1}.DMuc0;
        BCMObj.DLig0 = varargin{1}.DLig0;
        BCMObj.DMus0 = varargin{1}.DMus0;
        BCMObj.epsilon = varargin{1}.epsilon;
        BCMObj.Znodal = varargin{1}.Znodal;
        BCMObj.a_CT = varargin{1}.a_CT;
        BCMObj.a_TA = varargin{1}.a_TA;
        BCMObj.a_LC = varargin{1}.a_LC;
        BCMObj.mu = varargin{1}.mu;
        BCMObj.ml = varargin{1}.ml;
        BCMObj.mb = varargin{1}.mb;
        BCMObj.xu0 = varargin{1}.xu0;
        BCMObj.xl0 = varargin{1}.xl0;
        BCMObj.xb0 = varargin{1}.xb0;
        BCMObj.ku = varargin{1}.ku;
        BCMObj.kl = varargin{1}.kl;
        BCMObj.kc = varargin{1}.kc;
        BCMObj.kb = varargin{1}.kb;
        BCMObj.xi_01 = varargin{1}.xi_01;
        BCMObj.xi_02 = varargin{1}.xi_02;
        BCMObj.ParamSet_MuscleRules = varargin{1}.ParamSet_MuscleRules;
      else
        error('Incorrect input arguments! Valid options: ''male'', ''female'' or no argument at all!')
      end
      
%       if (strcmpi(BCMObj.sex,'male'))
%         BCMObj.Lg_init = BodyCoverModel.LREST_MALE;
%         BCMObj.Tg_init = BodyCoverModel.TREST_MALE;
%         BCMObj.DMuc_init = BodyCoverModel.DEPTH_MUC_MALE;
%         BCMObj.DLig_init = BodyCoverModel.DEPTH_LIG_MALE;
%         BCMObj.DMus_init = BodyCoverModel.DEPTH_MUS_MALE;
%       elseif (strcmpi(varargin{1},'female'))
%         BCMObj.Lg_init = BodyCoverModel.LREST_FEMALE;
%         BCMObj.Tg_init = BodyCoverModel.TREST_FEMALE;
%         BCMObj.DMuc_init = BodyCoverModel.DEPTH_MUC_FEMALE;
%         BCMObj.DLig_init = BodyCoverModel.DEPTH_LIG_FEMALE;
%         BCMObj.DMus_init = BodyCoverModel.DEPTH_MUS_FEMALE;
%       else
%         error('Acceptable ''Sex'' especification are ''male'' or ''female'' ')
%       end
      
      
    end
        
    function zetau = get.zetau(BCMobj) % CHECK!  : Equation 9
    % Function for computing zetau the damping factor for the upper mass 
      xu_f = BCMobj.xData(1);
      xu_col_f = BCMobj.xu_col;
      if (xu_f > xu_col_f)
        zetau = BCMobj.zetau0;
      else
        zetau = 2*BCMobj.zetau0 + 1.0; % 0.4;  
%         zetau = 2*zetau;
      end
    end
    
    function zetal = get.zetal(BCMobj) % CHECK! : Equation 10
    % Function for computing zetau the damping factor for the lower mass 
      xl_f = BCMobj.xData(2);
      xl_col_f = BCMobj.xl_col;
      if (xl_f > xl_col_f)
        zetal = BCMobj.zetal0;
      else
        zetal = 2*BCMobj.zetal0 + 1.0; % 0.4; 
%         zetal = 2*zetal;
      end
    end
    
    function hu_col = get.hu_col(BCMobj) % CHECK! Paragraph after Equation 6 b
    % Function for computing hu_col the linear spring constant during collision for the upper mass
      hu_col = 3*BCMobj.ku; % [N/m]
    end
    
    function hl_col = get.hl_col(BCMobj) % CHECK! Paragraph after Equation 6 b
    % Function for computing hl_col the linear spring constant during collision for the lower mass
      hl_col = 3*BCMobj.kl; % [N/m]
    end
    
    function du = get.du(BCMobj) % Equation 8a
    % Function for computing du the damping coefficient for the upper mass
      du = 2*BCMobj.zetau*sqrt(BCMobj.mu*BCMobj.ku); % [Ns/m]
    end
    
    function dl = get.dl(BCMobj) % Equation 8b
    % Function for computing dl the damping coefficient for the lower mass
      dl = 2*BCMobj.zetal*sqrt(BCMobj.ml*BCMobj.kl); % [Ns/m]
    end
    
    function db = get.db(BCMobj) % Equation 8c
    % Function for computing du the damping coefficient for the upper mass
      db = 2*BCMobj.zetab*sqrt(BCMobj.mb*BCMobj.kb); % [Ns/m]
    end
    
    function au = get.au(BCMobj) % Equation 11a
    % Function for computing au the glottal area for the upper portion
      deltau = (BCMobj.xData(1)-BCMobj.xu_col);
      au = max([0,2*BCMobj.Lg*deltau]); % [m^2]
    end
    
    function al = get.al(BCMobj) % Equation 11b
    % Function for computing al the glottal area for the lower portion
      deltal = (BCMobj.xData(2)-BCMobj.xl_col);
      al = max([0,2*BCMobj.Lg*deltal]); % [m^2]
    end
    
    function ag = get.ag(BCMobj) % ag from al and au
    % Function for computing ag the proyected glottal area
      ag = max([0, min([BCMobj.au, BCMobj.al])]); % [m^2]
    end
    
    function acont = get.acont(BCMobj) %HW check with Gabriel
    % Function for computing acont contact area
      acont = max([0, BCMobj.Lg*(BCMobj.Tu*(BCMobj.xData(1)<=0)+BCMobj.Tl*(BCMobj.xData(2)<=0))]); % [m^2]
    end
    
    function InitModel(BCMobj)
    % Function for initializing the dynamic state and simulation time index
    % prior to run the simulation of the Body Cover Model
      BCMobj.xData = zeros(9,1);
      BCMobj.xData(1) = BCMobj.xu0;
      BCMobj.xData(2) = BCMobj.xl0;
      BCMobj.xData(3) = BCMobj.xb0;
      
      BCMobj.n_IterCont = 0; % Simulation time index
    end
    
    function setState(BCMobj,Xstate)
    % Function for setting the dynamic state of the body cover model
      Xstate = Xstate(:);
      [row,col] = size(Xstate);
      if (row==9)&&(col==1)
        BCMobj.xData = Xstate;
      elseif (row==6)&&(col==1)
        BCMobj.xData = [Xstate; zeros(3,1)];
      elseif (row==3)&&(col==1)
        BCMobj.xData = [Xstate; zeros(6,1)];
      else
        error('Incorrect state vector dimension!')
      end
    end
    
    function setDrivingForceSolver(BCMobj,version)
    % Function for setting the rules version for computing the driving
    % forces due to aerodinamic pressure
      if ischar(version)&&strcmpi(version,'original')
        BCMobj.UseUpdated_ContactRules = false; 
        BCMobj.UseUpdated_AeroDrivingForces = false; 
      elseif ischar(version)&&strcmpi(version,'new')
        BCMobj.UseUpdated_ContactRules = true; 
        BCMobj.UseUpdated_AeroDrivingForces = true; 
      else
        error('Incorrect ''version'' string variable! Valid options are: ''original'', ''new''.')
      end
    end
    
    function setLgInit(BCMobj,LgVal)
    % Function for setting the rest length of the vocal fold
      BCMobj.Lg_init = LgVal;
      BCMobj.scalingVocalFold;
    end
    
    function setNonLinMode(BCMobj,ModeState)
    % Function for setting on (ModeState = true) or off (ModeState = false) 
    % the non-linear terms in the rules involved in the Body-Cover 
    % vocal fold model.
      if strcmp(ModeState,'on')
          BCMobj.NonLinMode = 1;
      elseif strcmp(ModeState,'off')
          BCMobj.NonLinMode = 0;
      else
          error('Available options: ''on'' for using non-linear terms, ''off'' for disregardin non linear terms.')
      end
    end
    
    function agcalc = calcGlottalArea(BCMobj)
    % Function for computing the glottal area
      agcalc = BCMobj.ag;
    end
    
    function aucalc = calcAreagu(BCMobj)
    % Function for computing the glottal area for the upper mass
      aucalc = BCMobj.au;
    end
    
    function alcalc = calcAreagl(BCMobj)
    % Function for computing the glottal area for the lower mass
      alcalc = BCMobj.al;
    end
    
    function acontcalc = calcContactArea(BCMobj)
    % Function for computing the glottal area
      acontcalc = BCMobj.acont;
    end
    
    % Functions defined on separate files : Equations of Motion
    [Fku,Fkl,Fkb,Fkc] = ElasticForces(BCMobj) % Equation 2/3/4
    
    [Fdu,Fdl,Fdb] = DampingForces(BCMobj) % Equation 7a/7b/7c
    
    [Fu_col,Fl_col] = CollisionForces(BCMobj) % Equation 6a/6b
    
    [Feu, Fel] = AeroPressure2DrivingForces(BCMobj, Ps, Pe, varargin) % Equation 23a/23b
    
    setSimulationParameter(BCMobj,Param)
    
    Simulate(BCMobj, Ps, Pe)
    
    scalingVocalFold(BCMobj,varargin)
        
  end
% - END CLASS DEFINITION -
end
