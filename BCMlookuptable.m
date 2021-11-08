%function BCMdata = BCM(a_CT,a_TA,Ps,priorF0,priorSPL)
clear all; close all; clc;
%% This is a simple example of how to use the BodyCoverClass

% Simulation parameters
fs = 50000; % [Hz] Sampling frequency HW leave at high resolution for the time being . However, we only need 11025 Hz.
%Ps = 800; % [Pa] Subglottal pressure
Pe = 0; % [Pa] Supraglottal pressure
T_Sim = 0.400; % [s] Simulation time  % 400 ms
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
%a_CT=0.25;
%a_TA=0.1;
a_LC=0.5;

% fundermenta frequency range
f_min = 10;
f_max = 900;

% initialize output arrays
output_param_Ug  = {};
output_param_Ag  = {};
output_param_x_aux  = {};

% nested for loop here
Ps = 10:100:2010; % pressure
a_TA = 0:0.02:1;
a_CT = 0:0.02:1;
timerset =tic;
for i = 9%1:1:length(Ps) % pressure
    for j = 7% 1:1:length(a_TA) % TA activation
        
        for k =7% 1:1:length(a_CT) % CT activation

          
          % Simulation CODE

          % Initialize
          VFObj.InitModel;
          x_out(:,1) = VFObj.xData;
          ag_out(1,1) = VFObj.ag;
          ag_out(2,1) = VFObj.au;
          ag_out(3,1) = VFObj.al;


          parfor simcont = 1: N_Sim
              MuscleActivation.Rule2BodyCoverParameters(VFObj,a_CT(k),a_TA(j),a_LC);
              VFObj.Simulate(Ps(i), Pe)
              x_aux(simcont,:) = VFObj.xData;
              Ag(simcont) = VFObj.ag;
              Qg(simcont) = getFlow(Ag(simcont),Ps(i),constFlow);
          end

          % Use the second part of Ug for F0 and SPL calculations
          Qg_stable =  Qg(end/2:end);
          F0 = measures_getf0(Qg_stable,fs,f_max,f_min);

          % calculate Pg from Qg_stable 
          impedence_const = 45000;%1000;
          Pg = impedence_const.*gradient(Qg_stable); 
          SPL = measures_getSPL(Pg,fs) ;

          % Save all information in matrix format % for ease of graphing

%           output_param(i,j,k).Ug= Qg';
%           output_param(i,j,k).Ag= Ag';
%           output_param(i,j,k).F0 = F0;
%           output_param(i,j,k).SPL = SPL;
%           
          
          
          
          output_param(j,k).Ug= Qg';
          output_param(j,k).Ag= Ag';
          output_param(j,k).F0 = F0
          output_param(j,k).SPL = SPL
          
          
         disp(['i=', num2str(i),' j=', num2str(j),' k=', num2str(k)]); toc(timerset);
        end
    end
    
end
%save('output_param_04212021_1.mat','output_param');
%end


