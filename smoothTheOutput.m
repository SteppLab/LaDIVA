load('PermutationsOut_19-Nov-2015_02-18-58 9hr_overnight.mat',...
    'input', 'output')

if 1 % scale data to n1-1 range
    % go from 0-1 first...
    input(:, 2) = (input(:, 2) - 2000) / 20000;
    input(:, 3) = input(:, 3) / .03;
    
    input = input * 2 - 1;
end

nInput = size(input, 1);

%{
nEpochs = 5;
smRad = .05;

tmpOut = output;
for n = 1:nEpochs
    fprintf('epoch %i of %i\n', n, nEpochs);
    newOut = zeros(size(tmpOut));
    
    for i = 1:size(input, 1)
        idx = randsample(nInput, round(nInput/2));
        subinput = input(idx, :);
        suboutput = output(idx, :);
        
        curr = input(i, :);
        d = sqrt(sum(bsxfun(@minus, subinput, curr).^2, 2));
        newOut(i, :) = mean(suboutput(d < smRad, :));
    end
    
    tmpOut = newOut;
    
end

%}

%input1 = reshape(input(:, 1), [61, 37, 51]);

output1 = reshape(output(:, 1), [61, 37, 51]);
output2 = reshape(output(:, 2), [61, 37, 51]);
output3 = reshape(output(:, 3), [61, 37, 51]);

output1 = smooth3(output1);
output2 = smooth3(output2);
output3 = smooth3(output3);

outputSm = [reshape(output1, [nInput, 1]), reshape(output2, [nInput, 1]),...
    reshape(output3, [nInput, 1])];

%% don't include a data point that is completely surround by F0 = 0
keep = false(size(output1));
for x = 1:61
    for y = 1:37
        for z = 1:51
            if output1(x, y, z) > 0
                keep(x, y, z) = true; continue; 
            end
            if x > 1 && output1(x-1, y, z) > 0
                 keep(x, y, z) = true; continue;
            end
            if x < 61 && output1(x+1, y, z) > 0
                keep(x, y, z) = true; continue; 
            end
            if y > 1 && output1(x, y-1, z) > 0
                keep(x, y, z) = true; continue; 
            end
            if y < 37 && output1(x, y+1, z) > 0
                keep(x, y, z) = true; continue; 
            end
            if z > 1 && output1(x, y, z-1) > 0
                keep(x, y, z) = true; continue;  
            end
            if z < 51 && output1(x, y, z+1) > 0
                keep(x, y, z) = true; continue;  
            end
        end
    end
end

keepC = reshape(keep, [nInput, 1]);
inputK = input(keepC, :);
outputK = outputSm(keepC, :);

%% now fit a GMM to the smoothed outputs and constrained inputs
Art = inputK;
voice = outputK;

% GMM fit
xy=Art;
sxy=std(xy);
xy=xy(:,sxy>0);
ncl=32;
maxiter=1000;
options=struct('Display','iter','MaxIter',maxiter);
gmfit_art=gmdistribution.fit(xy,ncl,'Start','randSample','Replicates',1,'SharedCov',true,'Regularize',1e-5,'Options',options);

% train the radial basis functions
x=gmfit_art.posterior(xy);
x=reshape(bsxfun(@times, x, permute([Art,ones(size(Art,1),1)],[1,3,2])),size(x,1),[]);
B_voice=pinv(x'*x)*x'*voice;
r_fmt=corrcoef([voice,x*B_voice]);

% keep the relevant statistics
iSigma=pinv(gmfit_art.Sigma);
MU=gmfit_art.mu;
P=gmfit_art.PComponents.';
fmfit=struct('mu',MU,'iSigma',iSigma,'p',P,'beta_voice',B_voice');

%% load saved fit
load('divaVoiceSmooth.mat', 'fmfit');
MU = fmfit.mu;
iSigma = fmfit.iSigma;
P = fmfit.p;
B_voice = fmfit.beta_voice';

%% how good is the gmm?
pd = pdf(gmfit_art, input);
pd3 = reshape(pd, [61, 37, 51]);

figure(400); clf;
subplot(1, 2, 1);
scatter3(inputK(:, 1), inputK(:, 2), inputK(:, 3));
xlabel('in1'); ylabel('in2'); zlabel('in3');

subplot(1, 2, 2);
patch(isocaps(pd3, .5), 'FaceColor', 'interp', 'EdgeColor', 'none');
xlabel('in1'); ylabel('in2'); zlabel('in3');
title('pdf');

%% test the fit
voice_fit=zeros(size(output));
for n1=1:size(input,1),
    dx=bsxfun(@minus,input(n1,:),MU);
    p=-sum((dx*iSigma).*dx,2)/2;
    p=P.*exp(p-max(p));
    p=p/sum(p);
    px=p*[input(n1,:),1];
    voice_fit(n1,:)=px(:)'*B_voice;
end

voice_fitK = voice_fit(keepC, :);
good = output(:, 1) > 0;
good = keepC;

%% plot too many things
figure(100); clf;
labels = {'CT tension', 'Lung pressure', 'VF distance'; 'F0', 'Intensity', 'HNR'};
for r = 1:3
    for c = 1:3
        subplot(3, 3, (r-1)*3 + c);
        plot(input(:, r), voice_fit(:, c), 'rx');
        xlabel(labels{1, r});
        ylabel(labels{2, c});
    end
end

%% how good is this fit?
% this evaluates the accuracy at a given tension
badInput = input(~keepC, :);
badOutput = output(~keepC, :);
badvoice_fit = voice_fit(~keepC, :);

ctTension = 1;
blah = inputK(:, 1) == ctTension; % which tension do you want to look at?
blah2 = badInput(:, 1) == ctTension;

figure(44); clf;
subplot(2, 2, 1);
hist(voice_fitK(blah, 1), 100);
title(sprintf('fit of good inputs: ctTension = %1.2f', ctTension));
subplot(2, 2, 3);
hist(outputK(blah, 1), 100);
title(sprintf('good outputs: ctTension = %1.2f', ctTension));
xlabel('F0');

subplot(2, 2, 2);
hist(badvoice_fit(blah2, 1), 100);
title(sprintf('bad inputs: ctTension = %1.2f', ctTension));
subplot(2, 2, 4);
hist(badOutput(blah2, 1), 100);
title(sprintf('fit of bad inputs: ctTension = %1.2f', ctTension));
xlabel('F0');

%% how good is this fit?
% this evaluates the accuracy at a given pressure
badInput = input(~keepC, :);
badOutput = output(~keepC, :);
badvoice_fit = voice_fit(~keepC, :);

pressure = .8;
blah = inputK(:, 2) == pressure; % which tension do you want to look at?
blah2 = badInput(:, 2) == pressure;

figure(45); clf;
subplot(2, 2, 1);
hist(voice_fitK(blah, 1), 100);
title(sprintf('fit of good inputs: pressure = %1.2f', pressure));
subplot(2, 2, 3);
hist(outputK(blah, 1), 100);
title(sprintf('good outputs: pressure = %1.2f', pressure));
xlabel('F0');

subplot(2, 2, 2);
hist(badvoice_fit(blah2, 1), 100);
title(sprintf('bad inputs: pressure = %1.2f', pressure));
subplot(2, 2, 4);
hist(badOutput(blah2, 1), 100);
title(sprintf('fit of bad inputs: pressure = %1.2f', pressure));
xlabel('F0');

%% how good is this fit?
% this evaluates the accuracy at a given distance
badInput = input(~keepC, :);
badOutput = output(~keepC, :);
badvoice_fit = voice_fit(~keepC, :);

distance = -.6;
blah = inputK(:, 3) == distance; % which tension do you want to look at?
blah2 = badInput(:, 3) == distance;

figure(46); clf;
subplot(2, 2, 1);
hist(voice_fitK(blah, 1), 100);
title(sprintf('fit of good inputs: distance = %1.2f', distance));
subplot(2, 2, 3);
hist(outputK(blah, 1), 100);
title(sprintf('good outputs: distance = %1.2f', distance));
xlabel('F0');

subplot(2, 2, 2);
hist(badvoice_fit(blah2, 1), 100);
title(sprintf('bad inputs: distance = %1.2f', distance));
subplot(2, 2, 4);
hist(badOutput(blah2, 1), 100);
title(sprintf('fit of bad inputs: distance = %1.2f', distance));
xlabel('F0');

%% plot the outputs
figure(77); clf;
subplot(3,1,1); hold on;
plot(input(good,1),output(good,1),'bo');
plot(input(~good,1),output(~good,1),'rx');
xlabel('CT tension'); ylabel('F0');
yl1 = ylim;
subplot(3,1,2); hold on;
plot(input(good,2),output(good,2),'bo');
plot(input(~good,2),output(~good,2),'rx');
xlabel('Lung pressure'); ylabel('Intensity');
yl2 = ylim;
subplot(3,1,3); hold on;
plot(input(good,3),output(good,3),'bo');
plot(input(~good,3),output(~good,3),'rx');
xlabel('VF distance'); ylabel('HNR');
yl3 = ylim;

figure(78); clf;
subplot(3,1,1); hold on;
plot(input(good,1),voice_fit(good,1),'bo');
plot(input(~good,1),voice_fit(~good,1),'rx');
xlabel('CT tension'); ylabel('F0');
ylim(yl1);
subplot(3,1,2); hold on;
plot(input(good,2),voice_fit(good,2),'bo');
plot(input(~good,2),voice_fit(~good,2),'rx');
xlabel('Lung pressure'); ylabel('Intensity');
ylim(yl2);
subplot(3,1,3); hold on;
plot(input(good,3),voice_fit(good,3),'bo');
plot(input(~good,3),voice_fit(~good,3),'rx');
xlabel('VF distance'); ylabel('HNR');
ylim(yl3);

%% how should I scale down these large negative values?
% diva only needs to know which way to move to get back to positive... it
% doesn't need to deal with -100000
F0 = voice_fit(voice_fit(:, 1) < 0, 1);
in = voice_fit(voice_fit(:, 2) < 0, 2);
hnr = voice_fit(voice_fit(:, 3) < 0, 3);

figure(72); clf;
subplot(3, 1, 1);
hist(F0/1500, 100);
subplot(3, 1, 2);
hist(-log(-in), 100);
subplot(3, 1, 3);
hist(hnr/100, 100);
