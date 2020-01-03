function [all_pulses, this_pulse, f0] = LeTalker_run_continuously(ct,ta,pressure,voicing, restartSimulation) %[p, r] = LeTalker_run_continuously(p,c,N,Fs,ctvect,tavect,plvect)
%LeTalker version 1.21
%LeTalker - means "Lumped element Talker" and is the control code
%for running the three-mass model of Story and Titze (1995) and Titze and Story (2002) with sub and supra glottal
% airways;
% Author: Brad Story, Univ. of Arizona
%Date: 06.12.2013

%maybe dumb, but basically leave everything in memory
persistent r p c n dt SUPRA SUB neq F xnois ynois Unois;
persistent    csnd rho rhoc mu  PI  bprev fprev Pox tme R;
persistent alpha  tcos fcos ug bigN Fs ieplx f b a jtrch itrch jmoth;
persistent  N last_n cb ca save_rescaled; 

if isempty(p) || restartSimulation
% % %     keep_track_tmp = 0;
    InitializeLeTalker
    p = rules_consconv(p,c);
    n=1; 
    bigN = 100000;
    update_fs = 200;
    Fs = 11025;
    N = ceil(Fs / update_fs);
    last_n = 0;
    %figure; 
    %load filter coefficients for band pass filtering glottal noise
    load filter_coeffs;

    %If SUB=1 the subglottal system is included; SUB=0 will take it out.
    SUB = 0; %p.SUB;
    if(SUB ==0)
        p.ar(end) = 1000000;
    end

    %If SUPRA=1 the supraglottal system is included; SUPRA=0 will take it out.
    SUPRA = 0; %p.SUPRA;
    if(SUPRA == 0)
        p.ar(1) = 1000000;
    end


    neq = 6;
    F = zeros(1,neq);
    r.theta = zeros(1,bigN);
    r.x = zeros(1,bigN);
    r.xb = zeros(1,bigN);

    p.tan_th0 = p.xc/p.T;
    p.th0 = atan(p.tan_th0);

    for i=1:6
        F(i) = 0;
    end;

    %For use with Titze 1 below-----
    F(1) = 0;
    F(2) = 0;
    F(3) = 0;
    %For use with Titze 2 below-----
    % F(1) = p.x01;
    % F(2) = p.x02;
    % F(3) = p.xb0;


    %-----
    a = p.ar(1,:);
    Nsect = length(a);
    ieplx = 1;
    jmoth = 44;
    itrch = jmoth+1;
    jtrch = Nsect;
    f = 	zeros(1,Nsect);
    b = 	zeros(1,Nsect);
    r1 = 	zeros(1,Nsect);
    r2 = 	zeros(1,Nsect);
    aprev = zeros(1,jmoth);
    atot = zeros(1,jmoth);

    r.x1 = []; %lower mass disp
    r.x2 = []; %upper mass disp
    r.xb = []; %body mass disp
    r.po = []; %output pressure
    r.ug = []; %glottal flow
    r.ga = []; %glottal area

    %---Initialize vectors for noise generation--
    xnois = zeros(1,5);
    ynois = zeros(1,5);
    Unois = 0;

    %vocal tract attenuation - same everywhere
    %alpha = .98*ones(1,Nsect);
    csnd =35000;
    rho  =0.00114;
    rhoc =rho*csnd;
    mu = 0.000186;
    PI = pi;
    %PI2 = 2*PI;

    bprev = 0;
    fprev = 0;
    Pox = 0.0;

    dt = 1/Fs;
    tme = 0.0;
    R = 128.0/(9.0*pi.^2);

    %--initialize----------
    p.a1 = 2*p.L*F(1);
    p.a2 = 2*p.L*F(2);
    p.ga = max(0,min(p.a1,p.a2));

    %vocal tract attenuation - based on xsect area
    alpha = 1-p.vtatten./sqrt(a);

    %posterior glottal gap
    pgap = 0;
    p.psg = p.ps;

    %------

    tcos = 0.002; %%DON'T RAMP ### MJC
    
    save_rescaled = [];
end

    %pressure is from [-1 1]; rescale to LeTalker range [2000 to 20,000] 
    %The general case (when you have a value c between a and b and you want a value x between y and z), x is calculated as follows:
    %x := (c - a) * (z - y) / (b - a) + y
    rescaledPressure = (pressure - 0) * (20000 - 2000) / (1 - 0)  + 2000; %a=-1, b=1; y = 2000 z=20000 
    %rescaledPressure = (pressure) * (20000 - 0) + 0; %##For now, min pressure = 7000; otherwise no sound
    rescaledVoicing = (voicing - 0) * (0 - .03) / (1 - 0) + .03;   %a=-1, b=1; y = .05 z=0 
    rescaledCT = (ct - 0) * (1 - 0) / (1 - 0) + 0;  %-1 to 1 vs 0 to 1 
    %disp(rescaledVoicing);
    p.Noise = 1; %%###%MJC tmp
    save_rescaled = [save_rescaled; [rescaledCT rescaledPressure rescaledVoicing]];
% % %     keep_track_tmp = keep_track_tmp + 1;
% % %     if ~mod(keep_track_tmp, 10); disp(rescaledVoicing); end;
for this_n=1:N
    n = last_n + this_n; %disp(n);
    t = (n-1)*dt;
    
    %-- check for time-varying muscle activation and pressure ramp
    p.ata = ta; 
%     if(isempty(tavect)==0)
%         p.ata=tavect(n);
%     end
    p.act = rescaledCT; 
% %     if(isempty(ctvect)==0)
% %         p.act=ctvect(n);
% %     end
%     if(isempty(plvect)==0)
        p.ps=rescaledPressure; %plvect(n);
        if(t <= tcos)                                 % ramp the driving pressure
            fcos = 1/(2*tcos);
            p.ps = p.ps*(cos(2*pi*fcos*t - pi)+1)/2;  % subglottal pressure
        end
        
        p.x02 = rescaledVoicing; %###.01; %voicing;
        p.x01 = rescaledVoicing; %### .011; %voicing; 

%     else
%         if(t <= tcos)                                 % ramp the driving pressure
%             fcos = 1/(2*tcos);
%             p.ps = p.PL*(cos(2*pi*fcos*t - pi)+1)/2;  % subglottal pressure
%         else
%             p.ps = p.PL;
%         end;
%     end
    
    %-------------------------
    
    p = rules_consconv(p,c);
    
    %This is internally necessary because p.ps was changed to p.psg in calc_pressures.m to accomodate the trachea wave
    %propagation for other implementations
    p.psg = p.ps;
  
    
    
    %================Titze 1 ==========================
    
    znot = p.zn/p.T;
    p.xn = (1-znot)*(p.x01+F(1))+znot*(p.x02 + F(2));
    p.tangent = (p.x01-p.x02 + 2*(F(1)-F(2)))/p.T;
    p.x1 = p.xn - ( -p.zn)*p.tangent;
    p.x2 = p.xn - (p.T - p.zn)*p.tangent;
    %====================================================
    
    
    %================Titze 2 ==========================
    %  p.tan_th0 = (p.x01-p.x02)/p.T;
    %  p.th0 = atan(p.tan_th0);
    %  p.tangent = tan(p.th0 + F(1));
    %  p.xn = p.x02 + (p.T - p.zn)*p.tan_th0 + F(2);
    %  p.x1 = p.xn - ( -p.zn)*p.tangent;
    %  p.x2 = p.xn - (p.T - p.zn)*p.tangent;
    %====================================================
    
    % Vocal fold area calculations
    p.zc = min(p.T, max(0., p.zn+p.xn/p.tangent));
    p.a1 = max(p.delta, 2*p.L*p.x1);
    p.a2 = max(p.delta, 2*p.L*p.x2 );
    p.an = max(p.delta, 2*p.L*p.xn);
    p.zd = min(p.T, max(0, -0.2*p.x1/p.tangent));
    p.ad = min(p.a2, 1.2*p.a1);
    
    % contact area calculations
    if(p.a1 > p.delta & p.a2 <= p.delta)
        p.ac = p.L*(p.T - p.zc);
    elseif(p.a1 <= p.delta & p.a2 > p.delta)
        p.ac = p.L*p.zc;
    elseif(p.a1 <= p.delta & p.a2 <= p.delta)
        p.ac = p.L*p.T;
    end;
    
    
    %Calculate driving pressures
    p.pe = f(ieplx)+b(ieplx);
    p = calc_pressures(p);
    
    %Integrate EOMs with Runge-Kutta
    F = rk3m(neq,F,tme,dt,p);
    
    %Glottal area calculation
    p.ga = max(0,min(p.a1,p.a2));
    
    % Calculate the glottal flow;
    ug = calcflow(p.ga,f(jtrch),b(ieplx),a,csnd,rhoc);
    
    
    %=========noise generator====================
    %Added on 10.04.12 - Reynolds number calculation and option for adding
    %noise to simulate turbulence-generated sound, based on Samlan and
    %Story (2011) JSLHR
    if(p.Noise == 1)
        
        RE2 = (ug * rho/ (mu*p.L))^2;
        RE2b = 1440000; %Corresponds to RE = 1200
        Anois = 0.000004*(RE2-RE2b);
        Unois = 0;
        
        tmp1 = cb(1)*xnois(5)+cb(2)*xnois(4) + cb(3)*xnois(3) + cb(4)*xnois(2) + cb(5)*xnois(1);
        tmp2 = ca(2)*ynois(4) + ca(3)*ynois(3) + ca(4)*ynois(2) + ca(5)*ynois(1);
        ynois(5) = tmp1-tmp2;
        
        if(Anois > 0)
            Unois = Anois * ynois(5);
        end
        
        ug = ug + Unois;
        
        for i=1:4
            xnois(i) = xnois(i+1);
            ynois(i) = ynois(i+1);
        end
        xnois(5) = rand(1)-0.5;
    end
    
    
    %   ============== wave propagation in vocal tract ======================= */
    
    D = a(1:jtrch-1) + a(2:jtrch);
    r1 = (a(1:jtrch-1) - a(2:jtrch))./D;
    r2 = -r1;
    
    
    % ---even sections in trachea--- */
    
    %f(ieplx) = alpha*b(ieplx) + p.u(n)*(rhoc/a(ieplx));
    f(ieplx) = alpha(ieplx)*b(ieplx) + ug*(rhoc/a(ieplx));
    
    if(SUPRA == 0)
        f(ieplx) = 0;
    end
    
    if(SUB == 1)
        p.psg = 2.*f(jtrch)*alpha(jtrch) - ug*(rhoc/a(jtrch));
        f(itrch)= 0.9*p.ps - 0.8*b(itrch)*alpha(itrch);
    else
        p.psg = p.ps;
        f(jtrch)= p.ps;
    end;
    
    f  = alpha.*f;
    b = alpha.*b;
    
    % ---even sections in trachea--- */
    if(SUB==1)
        
        Psi = f(itrch+1:2:jtrch-2).*r1(itrch+1:2:jtrch-2)  + b(itrch+2:2:jtrch-1).*r2(itrch+1:2:jtrch-2);
        b(itrch+1:2:jtrch-2) = b(itrch+2:2:jtrch-1) + Psi;
        f(itrch+2:2:jtrch-1) = f(itrch+1:2:jtrch-2) + Psi;
        b(jtrch) = f(jtrch)*alpha(jtrch) -ug*rhoc/a(jtrch);
    end
    
    % ---even sections in supraglottal--- */
    
    Psi = f(ieplx+1:2:jmoth-2).*r1(ieplx+1:2:jmoth-2)  + b(ieplx+2:2:jmoth-1).*r2(ieplx+1:2:jmoth-2);
    b(ieplx+1:2:jmoth-2) = b(ieplx+2:2:jmoth-1) + Psi;
    f(ieplx+2:2:jmoth-1) = f(ieplx+1:2:jmoth-2) + Psi;
    
    
    % odd sections */
    
    if(SUB ==1)
        Psi = f(itrch:2:jtrch-1).*r1(itrch:2:jtrch-1)  + b(itrch+1:2:jtrch).*r2(itrch:2:jtrch-1);
        b(itrch:2:jtrch-1) = b(itrch+1:2:jtrch) + Psi;
        f(itrch+1:2:jtrch) = f(itrch:2:jtrch-1) + Psi;
        
    end
    
    Psi = f(ieplx:2:jmoth-1).*r1(ieplx:2:jmoth-1)  + b(ieplx+1:2:jmoth).*r2(ieplx:2:jmoth-1);
    b(ieplx:2:jmoth-1) = b(ieplx+1:2:jmoth) + Psi;
    f(ieplx+1:2:jmoth) = f(ieplx:2:jmoth-1) + Psi;
    
    %------------- Lip Radiation -----------  */
    am = a(jmoth);
    L = (2.0/dt)*8.0*am/(3.0*PI*csnd);
    a2 = -R - L + R*L;
    a1 = -R + L - R*L;
    b2 = R + L + R*L;
    b1 = -R + L + R*L;
    b(jmoth) = (1/b2)*(f(jmoth)*a2+fprev*a1+bprev*b1);
    Pout = (1/b2)*(Pox*b1 + f(jmoth)*(b2+a2) + fprev*(a1-b1));
    Pox = Pout;
    
    bprev = b(jmoth);
    fprev = f(jmoth);
    
    % --Complete reflection at lips------
    % b(jmoth) = -f(jmoth);
    %Pout = b(jmoth)+f(jmoth);
    %---------------------------------------- */
    
    tme  = tme +dt;
    
    
    %Assign output signals to structure--------
    r.po(n) = Pout;
    r.ug(n) = ug;
    r.nois(n) = Unois;
    r.ga(n) = p.ga;
    r.ad(n) = p.ad;
    r.x1(n) = p.x1;
    r.x2(n) = p.x2;
    r.xb(n) = F(3);
    r.ps(n) = p.psg;
    r.pi(n) = p.pe;
    r.f1(n) = p.f1;
    r.f2(n) = p.f2;
    
    
    
    if(max(isnan(r.x1)) == 1 | max(isnan(r.x2)) | max(isnan(r.xb)) )
        r.x1 = zeros(bigN,1);
        r.x2 = zeros(bigN,1);
        break;
    elseif(max(isinf(r.x1)) == 1 | max(isinf(r.x2)) | max(isinf(r.xb)) )
        r.x1 = zeros(bigN,1);
        r.x2 = zeros(bigN,1);
        break;
    end;
    
    
%   n = n+1;
%   if ~mod(n,100)
%       figure; plot(r.ga);
%       track_pulses = r.ga;
%       break;
%   end
end;

%plot(1:last_n+this_n, r.ga);
all_pulses = r.ug;%/10^3;  %glottal flow (r.ga is glottal area)
this_pulse = r.ug(last_n+1:end);%/10^3; %newly generated pulse

 %The general case (when you have a value c between a and b and you want a value x between y and z), x is calculated as follows:
 %x := (c - a) * (z - y) / (b - a) + y
%rescaledPulses = (all_pulses - 0) * (2 - 0) / (4816 -0 ) + 0; %a=0, b=4816; y =0 z= 1   %max glottal pressure from permutations
rescaledPulses = (all_pulses - 0) * (3 - 0) / (2000 -0 ) + 0; %a=0, b=4816; y =0 z= 1   %max glottal pressure from permutations



%this_pulse = rescaledPulse;
all_pulses = rescaledPulses;
disp(max(all_pulses));
this_pulse = rescaledPulses(last_n+1:end);
% % % figure(110);
% % % subplot(3,1,1);
% % % plot(r.ug);
% % % subplot(3,1,2); 
% % % plot(rescaledPulses); 
% % % subplot(3,1,3)
% % % plot(r.ug,rescaledPulses,'.');


%whatttttt is this doing ###
last_n = last_n+this_n; %reset 
tmp = r.x1;
n1 = round(0.4*length(tmp));
n2 = length(tmp)-round(0.1*length(tmp));
tmp = tmp(n1:n2);
tmp = tmp;
if(min(tmp)>0);
    tmp = tmp-0.1*max(tmp);
end

[tmpf0,tme] = zerocross(tmp,Fs);

r.f0 = round(mean(tmpf0));
f0 = r.f0;
%figure(123); plot(save_rescaled);
%save_rescaled_tmp = save_rescaled;

