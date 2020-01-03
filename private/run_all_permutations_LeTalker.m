function run_all_permutations_LeTalker()

    InitializeLeTalker;
    fs = 11025;
    p.SUB = 0; %don't use letalker trachea
    p.SUPRA = 0; %don't use letalker vocal tract
    p.Noise = 1;
    N = 1000;
    
%% 22 hrs
% %     ct_steps = 0:.02:1;
% %     %ta_steps = 0:.01:1;
% %     pl_steps = 2000:200:20000;
% %     vfdist_steps = 0:.003:.18;
% %     time_est = length(ct_steps) * length(pl_steps) * length(vfdist_steps) * .28 / 60 / 60;

%% % %         3.44 hrs
% %     ct_steps = 0:.04:1;
% %     %ta_steps = 0:.01:1;
% %     pl_steps = 2000:400:20000;
% %     vfdist_steps = 0:.005:.18;
% %     time_est = length(ct_steps) * length(pl_steps) * length(vfdist_steps) * .28 / 60 / 60;

%% 30 mins
% %     ct_steps = 0:.05:1;
% %     %ta_steps = 0:.01:1;
% %     pl_steps = 2000:1000:20000;
% %     vfdist_steps = 0:.01:.18;
% %     time_est = length(ct_steps) * length(pl_steps) * length(vfdist_steps) * .28 / 60 / 60;

%% 11 mins
% % %     ct_steps = 0:.06:1;
% % %     %ta_steps = 0:.01:1;
% % %     pl_steps = 2000:2000:20000;
% % %     vfdist_steps = 0:.015:.18;
% % %     time_est = length(ct_steps) * length(pl_steps) * length(vfdist_steps) * .28 / 60 / 60;

    ct_steps = 0:.02:1;
    %ta_steps = 0:.01:1;
    pl_steps = 2000:500:20000;
    vfdist_steps = 0:.0005:.03;
    time_est = length(ct_steps) * length(pl_steps) * length(vfdist_steps) * .28 / 60 / 60;
    
%%
try
    tic
    %input = [tension, pressure, distance] 
    %output = [F0, intensity, noise to harmonic ratio?] high NHR = /h/ low = /u/
    idx = 1;
    input = nan(length(ct_steps)*length(pl_steps)*length(vfdist_steps),3);
    output = nan(length(ct_steps)*length(pl_steps)*length(vfdist_steps),3);
    glottal_pressure = nan(length(ct_steps)*length(pl_steps)*length(vfdist_steps),N);
    %glottal_noise = nan(length(ct_steps)*length(pl_steps)*length(vfdist_steps),N);
    %figure(7); hold on;
    for cidx = 1:length(ct_steps)
        %for tidx = 1:length(ta_steps)
            for pidx = 1:length(pl_steps)
                for vidx = 1:length(vfdist_steps)
                    %tic;
                    
                    if ~mod(idx,100)
                        disp(idx);
                    end
                    
                    ctvect = ones(1,N) * ct_steps(cidx);
                    tavect = ones(1,N) * .2;%ta_steps(tidx);
                    plvect = ones(1,N) * pl_steps(pidx);
                    p.x01 = vfdist_steps(vidx);
                    p.x02 = vfdist_steps(vidx);
                    [p, r] = LeTalker(p,c,N,fs,ctvect,tavect,plvect); 
                    %disp(toc);
                    
                    input(idx,:) = [ct_steps(cidx) pl_steps(pidx) vfdist_steps(vidx)];
                    output(idx,1) = r.f0; %F0
                    output(idx,2) = rms(r.ug(500:end)); %intensity - units??
                    [HNRlin, HNRdb] = calcNHR(r.ug, fs, r.f0);
                    %plot(NHR,HNR,'ko');
                    output(idx,3) =HNRlin; %NHR r.F0;
                    %output(idx,4) =HNRdb; %NHR r.F0;
                    glottal_pressure(idx,:) = r.ug;
                    %glottal_noise(idx,:) = r.nois;
                    idx = idx+1;
                end
            end
        %end
    end
    toc
    save(['PermutationsOut_' datestr(now, 'dd-mmm-yyyy_HH-MM-SS')]);
catch
    save(['PermutationsOut_' datestr(now, 'dd-mmm-yyyy_HH-MM-SS')]);
end



%93574589.08 s = 1559576.48466 mins = 25992.94141 hrs = 1083 days

figure(77)
subplot(3,1,1)
plot(input(:,1),output(:,1),'o');
xlabel('CT tension'); ylabel('F0');
subplot(3,1,2)
plot(input(:,2),output(:,2),'o')
xlabel('Lung pressure'); ylabel('Intensity');
subplot(3,1,3)
plot(input(:,3),output(:,3),'o')
xlabel('VF distance'); ylabel('HNR');


% figure(88)
% subplot(3,3,1)
% plot3(input(:,1),output(:,1),'o');
% xlabel('CT tension'); ylabel('F0');
% subplot(3,1,2)
% plot(input(:,2),output(:,2),'o')
% xlabel('Lung pressure'); ylabel('Intensity');
% subplot(3,1,3)
% plot(input(:,3),output(:,3),'o')
% xlabel('VF distance'); ylabel('HNR');



end
