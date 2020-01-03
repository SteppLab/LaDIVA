function [HNRlin, HNRdb] = calcNHR(y_raw, fs, f0)
    y_rms  = rms(y_raw);
    y_norm = y_raw/y_rms;
    
    if f0 == 0
        Mmax = 550;
        Mmin = 75;
    else
        Mmax = round(fs/(f0*.4)); %fs/sample = f0
        Mmin = round(fs/(f0*1.6));
    end
% % %     
% % % % %     toggle = 1; %If unknown f0, use hard coded values. Eventually make this a varagin
% % %     if toggle == 1 
% % % 
% % %     elseif toggle == 2
% % % 
% % %     end
% % %     
    y = y_norm - mean(y_norm);
    M = length(y);
    [correls, lags] = xcorr(y, 'coeff');
    correls(M) = 0;
    correls = correls(M+Mmin:M+Mmax); 
    
    Rmaxlag = max(correls);
    HNRlin = (Rmaxlag)/(1-Rmaxlag);
    HNRdb = 10*log10((Rmaxlag)/(1-Rmaxlag));
end