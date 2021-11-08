function [BCMdata]=BCMsynthesizer(a_CT,a_TA,Ps)
%% Sound Synthesizer with Body Cover Model
% Source signal : from BCM
% Filter from DIVA model : Digital filter 

% Set sampling rates %
fs = 50000;
art_fs = 200;

%resize Ps
maxPs = 2010; minPs=10;
Ps = Ps*(2010-10)+10; %resized up


% Make sure values of activation and PS are within required limits of BCM
a_CT = min(1,max(0,a_CT));
a_TA = min(1,max(0,a_TA));
Ps = min(2010,max(0,Ps)); %modified , was backwards %hasini 07/01/2020
% Manually setting values for BCM inputs %
% % Set 1
% a_CT = linspace(0.1,0.7,101);
% a_TA = linspace(0.1,0.7,101);
% Ps = 800*linspace(0.5,2,101);
% Set 2
%a_CT = 0.1*ones(101,1);
%a_TA = 0.1*ones(101,1);
%Ps = 800*ones(101,1);

% Setup BCM initial parameters %
Pe = 0; % [Pa] Supraglottal pressure
T_Sim = length(a_CT)/art_fs; % [s] Simulation time  % 400 ms (discarding initial 200ms due to trancient period)
N_Sim = floor(T_Sim*fs); % Simulation steps
Art_length = 0:1/art_fs:T_Sim-1/art_fs;
Sim_length = 0:1/fs:T_Sim-1/fs;

% interpolate 200Hz inputs to 500000Hz %
a_CTintp = interp1(Art_length,a_CT,Sim_length,'linear')'; a_CTintp = max(0,a_CTintp);
a_TAintp = interp1(Art_length,a_TA,Sim_length,'linear')'; a_TAintp = max(0,a_TAintp);
Psintp   = interp1(Art_length,Ps,Sim_length,'linear')'; Psintp = max(0,Psintp);

% Initialize Constant Flow parameters from Titze 1984 %
constFlow.PGO = 0 ; % Ap = 0 
constFlow.Ae = 3.0000e-04;     % Ae = 3 cm^2 = 3/10000 m^2
constFlow.As = 3.0000e-04;     % As = 3 cm^2 = 3/10000 m^2
constFlow.mu = 1.8600e-05;     % mu = 0.000186 g/(cm s) = 0.000186*(100/1000) kg/(m s)
constFlow.rho =  1.1400;       % rho = 0.00114 g/cm^3 = 0.00114*(1000000/1000) kg/m^3
constFlow.c = 350;             % c = 35000 cm/s = 350 m/s
constFlow.L = 0.0160;          % L = 1.6 cm = 1.6/100 m
constFlow.T = 3.5000e-05;      % T = 0.35 mm = 0.35/10000 m 
constFlow.solver = 'TITZE84M';

% Model calling and initialization %
VFObj=BodyCoverModel('male');
VFObj.setDrivingForceSolver('new')
VFObj.setSimulationParameter(fs)

% Muscle activation levels %
%a_CT=0.1;
%a_TA=0.1; 
a_LC=0.5;

% fundermenta frequency range %
f_min = 10;
f_max = 900;

% Simulation CODE %

  % Initialize
  VFObj.InitModel;
  x_out(:,1) = VFObj.xData;
  ag_out(1,1) = VFObj.ag;
  ag_out(2,1) = VFObj.au;
  ag_out(3,1) = VFObj.al;


  for simcont = 1: N_Sim
      MuscleActivation.Rule2BodyCoverParameters(VFObj,a_CTintp(simcont),a_TAintp(simcont),a_LC);
      VFObj.Simulate(Psintp(simcont), Pe)
      x_aux(simcont,:) = VFObj.xData;
      Ag(simcont) = VFObj.ag;
      Qg(simcont) = getFlow(Ag(simcont),Psintp(simcont),constFlow);
  end
  
  BCMdata.Ug = Qg';
  BCMdata.Ag = Ag';
  BCMdata.all_pulses = gradient(BCMdata.Ug); % Generate first derivative of glottal source to pass through filter

end
