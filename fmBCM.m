function [F0,SPL]=fmBCM(a_CT,a_TA,Ps,lookuptable1)
%% BCM forward model (single 1x1 value for each term above)
global priorCTindex

%resize Ps
maxPs = 2010; minPs=10;
Ps = Ps*(2010-10)+10; %resized up


% load lookup table 
%load('BCMlookuptable.mat');
[nPs,naTA,naCT] = size(lookuptable1.F0); % for now use the low res version as the interpolation needs work
aCTmin = 0; 
aTAmin = 0;
Psmin = lookuptable1.Psrange(1,1);

% for high resolution
aCTres = lookuptable1.aCTrange(1,2)-lookuptable1.aCTrange(1,1);
aTAres = lookuptable1.aTArange(1,2)- lookuptable1.aTArange(1,1);
Psres =  lookuptable1.Psrange(1,2)-lookuptable1.Psrange(1,1);

% for low resolution
% aCTres = 0.1;
% aTAres = 0.1;
% Psres = 50;

a_CT_round = aCTres.*round(a_CT./aCTres,0,'decimal');
a_TA_round = aTAres.*round(a_TA./aTAres,0,'decimal');
Ps_round = Psres.*round(Ps./Psres,0,'decimal');

% map inputs to outputs
CTindex = roundn(1+((a_CT_round-aCTmin)./aCTres),0);
TAindex = roundn(1+((a_TA_round-aTAmin)./aTAres),0);
Psindex = 1+((Ps_round-Psmin)./Psres);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets very coarse upper and lower bounds to the corrective motor command (feedforward signal)
if Psindex < 1,   Psindex = 1; end  
if Psindex > nPs,   Psindex = nPs; end

if CTindex < 1,   CTindex = 1; end
if CTindex > naCT,   CTindex = naCT; end

if TAindex < 1,   TAindex = 1; end
if TAindex > naTA,   TAindex = naTA; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(Psindex,2)> 1
    for lcount = 1:1:length(Psindex)
    F0(lcount,1) = lookuptable1.F0(Psindex(lcount),TAindex(lcount),CTindex(lcount));
    SPL(lcount,1) = lookuptable1.SPL(Psindex(lcount),TAindex(lcount),CTindex(lcount));
    end
else
    F0 = lookuptable1.F0(Psindex,TAindex,CTindex);
    
    SPL = lookuptable1.SPL(Psindex,TAindex,CTindex);
end
priorCTindex = CTindex ;
end