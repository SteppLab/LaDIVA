function SPL = measures_getSPL(signal,fs) % Input signal in Pa.
%   N = floor(length(signal)/windows);
%   m = zeros (1,N);
%   for count = 1:N
%     pos_from = (count-1)*windows+1;
%     pos_to = count*windows;
      m = 20*log10(rms(signal)/20e-6);
      %m = 20*log10(rms(signal)/0.5e-9);
%  end
  SPL = m;
end
