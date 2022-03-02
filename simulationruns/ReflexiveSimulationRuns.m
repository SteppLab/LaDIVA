%% changing reflexive gain levels and ploting results 

% trial_1 contains each run auditory outputs

 baselineHz = 133.9840;% Hz % Case B AA sustained vowel production vocal fo
% DIVA_x.pertresults.baselineHz = from baseline simulation
% DIVA_x.pertresults.perturbationHz
% -> time series of perturbation starting 550ms from sim start 
 DIVA_x.pertresults.perturbationHz  = 7.967;%Hz --> 99.9987 cents
% calcualtion: perturbation (cent) = 1200*log2((baseline_Hz+7.967)/baseline_Hz))
%
% DIVA_x.pertresults.baselineHz
% -> timeseries of an unperturbed trial auditory output vocal for model
% % DIVA_x.pertresults.meanreflex_f0HZ
% -> mean baseline response magnitude (calculation window same as reflexive)

% DIVA_x.pertresults.meanreflex_f0HZ
% -> mean reflexive analysis window vocal fo average
% -> analysis window = 120ms to 240ms from pert onset
% -> pert onset at 550ms from sim onset 550+(120:240) 
% -> 135:159 sample points at 200Hz (50ms windows)
DIVA_x.pertresults.Audfbgain_simsHz=[];
DIVA_x.pertresults.Audfbgain_sims_normdev_Hz=[];
DIVA_x.pertresults.Audfbgain_sims_normdev_cents=[];
DIVA_x.pertresults.Audfbgain_sims_time=DIVA_x.logs.time*1000 - 500;
DIVA_x.pertresults.perturbationHz=zeros(621,1);
DIVA_x.pertresults.perturbationHz(101:end,1) =perturbationHz.*ones(621-100,1);
%% consolidated results
% DIVA_x.pertresults.Audfbgain_sims_time  --> time from pertrubation onset (ms)
  sim_output= DIVA_x.pertresults.trial_1(:,2); % 2nd column is vocal fo
% for each sim
  audfbgain = 0.1;  %change value here after every run
% DIVA_x.pertresults.Audfbgain_simsHz --> all gains auditory output in Hz
  DIVA_x.pertresults.Audfbgain_simsHz = [DIVA_x.pertresults.Audfbgain_simsHz ,sim_output];
% DIVA_x.pertresults.Audfbgain_simslabels --> gain values for each column
  DIVA_x.pertresults.Audfbgain_simslabels = {{'audgain0.1'},{'audgain0.2'},{'audgain0.3'},{'audgain0.4'},{'audgain0.5'},{'audgain0.6'},{'audgain0.7'},{'audgain0.8'},{'audgain0.9'},{'audgain1'}};
% DIVA_x.pertresults.Audfbgain_sims_normdev_Hz --> auditory output - perturbation - baseline (Hz)
  response_Hz = sim_output - DIVA_x.pertresults.perturbationHzseries - DIVA_x.pertresults.baselineHzseries;
  DIVA_x.pertresults.Audfbgain_sims_normdev_Hz =[DIVA_x.pertresults.Audfbgain_sims_normdev_Hz, response_Hz];
% DIVA_x.pertresults.Audfbgain_sims_normdev_cents --> converted to cents
% calcualtion: response (cent) = 1200*log2((baseline_Hz+response(Hz))/baseline_Hz))
  response = 1200*log2((baselineHz+response_Hz)/baselineHz);
  DIVA_x.pertresults.Audfbgain_sims_normdev_cents=[DIVA_x.pertresults.Audfbgain_sims_normdev_cents,response];