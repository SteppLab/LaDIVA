function BCMdata = BCM(a_CT,a_TA,Ps)

%% This is how to use the BodyCoverClass
%resize Ps
maxPs = 2010; minPs=10;
%Ps = Ps*(2010-10)+10; %resized up

% Simulation parameters
fs = 50000; % [Hz] Sampling frequency HW leave at high resolution for the time being . However, we only need 11025 Hz.
%Ps = 800; % [Pa] Subglottal pressure
Pe = 0; % [Pa] Supraglottal pressure
T_Sim = 0.400; % [s] Simulation time  % 50 ms % input is given every 5ms - consider a 50ms window each time.
N_Sim = floor(T_Sim*fs); % Simulation steps

% Initialize Constant Flow parameters from Titze 1984 
constFlow.PGO = 0 ; % Ap = 0 
constFlow.Ae = 3.0000e-04;     % Ae = 3 cm^2 = 3/10000 m^2
constFlow.As = 3.0000e-04;     % As = 3 cm^2 = 3/10000 m^2
constFlow.mu = 1.8600e-05;     % mu = 0.000186 g/(cm s) = 0.000186*(100/1000) kg/(m s)
constFlow.rho =  1.1400;       % rho = 0.00114 g/cm^3 = 0.00114*(1000000/1000) kg/m^3
constFlow.c = 350;             % c = 35000 cm/s = 350 m/s
constFlow.L = 0.0160;          % L = 1.6 cm = 1.6/100 m
constFlow.T = 3.5000e-05;      % T = 0.35 mm = 0.35/10000 m 
constFlow.solver = 'TITZE84M';
%constFlow.delta_p = Ps; % 800 Pa Subglottal Pressure with a source only model

% Model calling and initialization
VFObj=BodyCoverModel('male');
VFObj.setDrivingForceSolver('new')
VFObj.setSimulationParameter(fs)

% Muscle activation levels
%a_CT=0.1;
%a_TA=0.1;
a_LC=0.5;

% fundermenta frequency range
f_min = 10;
f_max = 900;

% Simulation CODE

  % Initialize
  VFObj.InitModel;
  x_out(:,1) = VFObj.xData;
  ag_out(1,1) = VFObj.ag;
  ag_out(2,1) = VFObj.au;
  ag_out(3,1) = VFObj.al;


  for simcont = 1: N_Sim
      MuscleActivation.Rule2BodyCoverParameters(VFObj,a_CT,a_TA,a_LC);
      VFObj.Simulate(Ps, Pe)
      x_aux(simcont,:) = VFObj.xData;
      Ag(simcont) = VFObj.ag;
      Qg(simcont) = getFlow(Ag(simcont),Ps,constFlow);
  end
  
   Qg_stable =  Qg(end/2:end);
   F0 = measures_getf0(Qg_stable,fs,f_max,f_min);
%   if isnan(F0)
%        F0 = priorF0;
%   end
   % calculate Pg from Qg 
   impedence_const = 45000;
   Pg = impedence_const.*gradient(Qg_stable); 
   SPL = measures_getSPL(Pg,fs) ;
%   if isnan(SPL)
%        SPL = priorSPL;
%   end
%   
   % convert pressure to a value between 1 and 0
%   if SPL < 0
%      SPL = 1+ SPL./10000; 
%   end
%   
  BCMdata.F0 = F0
  BCMdata.SPL = SPL
  BCMdata.Ug = Qg';
  BCMdata.Ag = Ag';
  BCMdata.xData = x_aux;
 
 time=(0:1:length(BCMdata.Ag)-1).*(1000/fs);
%%
 figure; 
 subplot(2,1,1),plot(time,BCMdata.Ag); title ('BCM Output 1 = Ag'); xlabel('time (ms)'); ylabel('m^2');
 subplot(2,1,2), plot(time,BCMdata.Ug); title ('BCM Output 2 = Ug'); xlabel('time (ms)'); ylabel('m^3/s');
% 
 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');

end


