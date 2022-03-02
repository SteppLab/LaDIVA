function F0 = measures_getf0(signal,fs,f_max,f_min) % Input signal is glottal airflow 
  
      lim_inf = ceil(fs/(f_max));
      lim_sup = floor(fs/(f_min));
      U = xcov(signal','unbias'); % scales the raw covariance by 1/(M-abs(k)), where k is the index into the result.

      U = U(ceil(length(U)/2):end);
      %U = U(ceil(end/2):end);
      U = (U(lim_inf:lim_sup)-min(U(lim_inf:lim_sup)))/(max(U(lim_inf:lim_sup)) - min(U(lim_inf:lim_sup)));
      [M,P] = findpeaks(U); % M = peaks values, P =locations of peaks


      if isempty(P)
        F0 = NaN;
      else
        P = P(find(M>=0.9,1,'first'));
        if isempty(P)
          F0 = NaN;
        else
          F0 = fs/(P + lim_inf);
        end

        NFFT = pow2(nextpow2(length(signal)/4));
        [Pxx,Fxx] = pwelch(signal,NFFT,[],[],fs,'onesided');


        if ~isnan(F0)
          H = Pxx(find(Fxx>=F0,1,'first'));
          if (10*log10(max(Pxx)/H) > 80)
            F0 = NaN;
          end
        end
      end
   
  
end
