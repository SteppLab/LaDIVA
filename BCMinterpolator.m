%% BCM interpolator
%function [lookuptable] = BCMinterpolator(output_param)


[nPs,naTA,naCT] = size(output_param);
f0_calc = [output_param.F0];
f0_calc = reshape(f0_calc,[nPs,naTA,naCT]);% [Hz] fundamental frequency
SPL_calc = [output_param.SPL];
SPL_calc = reshape(SPL_calc,[nPs,naTA,naCT]);% [dB SPL] SPL

% f0_calc_1= f0_calc(1:11,:,:);
% SPL_calc_1= SPL_calc(1:11,:,:);
%f0_calc_2= f0_calc(11:21,:,:);
%SPL_calc_2= SPL_calc(11:21,:,:);

aCTmin = 0; aCTmax = 1;
aTAmin = 0; aTAmax = 1;
Psmin = 10; Psmax = 2010;%2010
%Psmin_2 = 1010;Psmax_2 = 2010;

% Original Resolutions
a_CTres_orig = 0.02; % a_CT : round to nearest 0.01
a_TAres_orig = 0.02; % a_TA : round to nearest 0.01
Psres_orig   = 100;   % Ps : round to nearest 10 Pa

a_CT_length = aCTmin:a_CTres_orig:aCTmax; %Z
a_TA_length = aTAmin:a_TAres_orig:aTAmax; %X
Ps_length   = Psmin:Psres_orig:Psmax;     %Y
%Ps_length_2   = Psmin_2:Psres_orig:Psmax_2;     %Y

% %Z = a_CT_length ; %Z
% %Y = a_TA_length ; %X
% %X = Ps_length;    %Y

 [X,Y,Z] = meshgrid(a_TA_length,Ps_length,a_CT_length);
%[X_2,Y_2,Z_2] = meshgrid(a_TA_length,Ps_length_2,a_CT_length);

% Desired lookup table Resolution
a_CTres = 0.01; % a_CT : round to nearest 0.01
a_TAres = 0.01; % a_TA : round to nearest 0.01
Psres   = 20;   % Ps : round to nearest 5 Pa

a_CT_length_new = aCTmin:a_CTres:aCTmax;
a_TA_length_new = aTAmin:a_TAres:aTAmax;
 Ps_length_new   = Psmin:Psres:Psmax;
%Ps_length_new_2   = Psmin_2:Psres:Psmax_2;

 [Xq,Yq,Zq] = meshgrid(a_TA_length_new,Ps_length_new,a_CT_length_new);
% [Xq_2,Yq_2,Zq_2] = meshgrid(a_TA_length_new,Ps_length_new_2,a_CT_length_new);

% 3D interpolation
 F0 = interp3(X,Y,Z,f0_calc,Xq,Yq,Zq,'cubic'); 
 SPL =interp3(X,Y,Z,SPL_calc,Xq,Yq,Zq,'cubic'); 
% F0_2 = interp3(X_2,Y_2,Z_2,f0_calc_2,Xq_2,Yq_2,Zq_2,'cubic');
% SPL_2 = interp3(X_2,Y_2,Z_2,SPL_calc_2,Xq_2,Yq_2,Zq_2,'cubic'); 

%F0(1:501,:,:) = F0_1;
% F0 = F0_2;
%SPL(1:501,:,:) = SPL_1;
% SPL = SPL_2;

lookuptable.F0 = F0;
lookuptable.SPL =SPL;

lookuptable.aCTrange = a_CT_length_new;
lookuptable.aTArange =  a_TA_length_new;
lookuptable.Psrange =Ps_length_new;
lookuptable.aCTrangeold = a_CT_length;
lookuptable.aTArangeold =  a_TA_length;
lookuptable.Psrangeold =Ps_length;


%end