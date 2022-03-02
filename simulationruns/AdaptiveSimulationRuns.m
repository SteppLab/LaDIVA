%% response transformation to cents

 DIVA_x.pertmagnitudeHz =  7.9670; % Case B
 baseline_Hz = 133.9840; % Case B
 pertdiffcent = 1200*log2((baseline_Hz + DIVA_x.pertmagnitudeHz)/baseline_Hz) ; % for pitch reflex
 DIVA_x.pertvalue = [zeros(1,DIVA_x.adaptphases(1)), linspace(0,DIVA_x.pertmagnitudeHz,DIVA_x.adaptphases(2)), DIVA_x.pertmagnitudeHz.*ones(1,DIVA_x.adaptphases(3)),zeros(1,DIVA_x.adaptphases(4))]; % HRW 09232020 added for pitch adaptation
 baselinetrial = DIVA_x.pertresults.trial_1;
adaptresponse=zeros(621,108);
 for i = 1:108
adaptresponse(:,i) =  eval(['DIVA_x.pertresults.trial_',num2str(i),'(:,1)-baselinetrial-DIVA_x.pertvalue(1,i).*ones(621,1);' ]);
 end

 DIVA_x.pertresults.adaptresponse=adaptresponse;
meanadapt_f0HZ_real = mean(DIVA_x.pertresults.adaptresponse(8:24,:),1); %40 ms to 120 ms
 meanadapt_f0cent_real=1200*log2((baseline_Hz+meanadapt_f0HZ_real)/baseline_Hz);
 
  DIVA_x.pertresults.meanadapt_f0HZ_real = [ DIVA_x.pertresults.meanadapt_f0HZ_real, meanadapt_f0HZ_real];
  DIVA_x.pertresults.meanadapt_f0cent_real =[DIVA_x.pertresults.meanadapt_f0cent_real,meanadapt_f0cent_real];