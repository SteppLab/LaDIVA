rescaledPulses = (all_pulses - 0) * (2 - 0) / (4816 -0 ) + 0; %a=0, b=4816; y =0 z= 1   %max glottal pressure from permutations
% this_pulse = rescaledPulse;
all_pulses = rescaledPulses;
this_pulse = rescaledPulses(last_n+1:end);
% % % figure(110);



    %The general case (when you have a value c between a and b and you want a value x between y and z), x is calculated as follows:
    %x := (c - a) * (z - y) / (b - a) + y
 c = 7840;
 a = 2000;
 b = 20000;
 y = 0;
 z = 1;
 x = (c - a) * (z - y) / (b - a) + y
 % 7840 pl is -.3511 in DIVA
 
 
 c = .2;
 a = 0;
 b = 1;
 y = 0;
 z = 1;
 x = (c - a) * (z - y) / (b - a) + y
 
 %.2 ct is  -.6  in DIVA
 
 
  
 c = .01;
 a = .05;
 b = 0;
 y = 0;
 z = 1;
 x = (c - a) * (z - y) / (b - a) + y
 %.01 dist = .6 voicing in DIVA
 
 
 
 load('PermutationsOut_19-Nov-2015_02-18-58 9hr_overnight.mat')
  set(0,'DefaultAxesFontSize',12)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % RESCALE DIVA _OUTPUTS_
 
 c = 7840;
 a = 2000;
 b = 20000;
 y = 0;
 z = 1;
 x = (c - a) * (z - y) / (b - a) + y
 % 7840 pl is -.3511 in DIVA
 
 
 c = .2;
 a = 0;
 b = 1;
 y = 0;
 z = 1;
 x = (c - a) * (z - y) / (b - a) + y
 
 %.2 ct is  -.6  in DIVA
 
 
  
 c = .01;
 a = .05;
 b = 0;
 y = 0;
 z = 1;
 x = (c - a) * (z - y) / (b - a) + y
 %.01 dist = .6 voicing in DIVA