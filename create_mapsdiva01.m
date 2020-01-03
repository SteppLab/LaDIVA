%STEPS=[1,2];
%% JI: this section just plots outputs from gridded inputs
if 0; % JI: this was originally 0
    idx_art=1:10;%1:5;%1:3;
    x=linspace(-1,1,32);
    for n=1:numel(idx_art)
        for m=1:numel(x)
            X=zeros(numel(idx_art),1);
            X(n)=x(m);
            [Aud(:,m,n),Som,Outline(:,m,n),af,filt(:,m,n)]=diva_synth([X;zeros(3,1)],'explicit');
        end
    end
    for k=idx_art
        subplot(2,5,k);plot(Aud(:,:,k)');
    end
    k=2;
    % JI: Outline plots the profile of the face/ vocal tract and filter
    subplot(121);plot(Outline(:,:,k));set(gca,'xdir','reverse'); axis equal;
    subplot(121);hold on; plot(Outline(:,end,k),'k','linewidth',2);hold off
    subplot(122);plot(Aud(:,:,k)');
    %subplot(122);plot(abs(filt(:,:,k)));
    %subplot(122);hold on; plot(abs(filt(:,end,k)),'k','linewidth',2);hold off
end

%% JI: this section runs diva_synth on random inputs and saves 
% note1: for a new vt model, save variable vt in rtart.mat to diva_synth.mat
% note2: remember to clear diva_synth.m if you have changed the
% diva_synth.mat file before running this script
if 1; %JI: originally - any(STEPS==1),
    %note: set reflection coefficients to .99 in diva_synth to produce narrower formants
    %vt.Average(353:354)=nan;
    min_af=.05;
    idx_art=1:10;%1:5;%1:3;
    nsamples=1e5;
    %sampling=12;%64;
    scale=1;%.5;%1.25;
    numf=3;
    samplesperperiod=1000;
    fs=11025;
    numdirs=numel(idx_art);
    X=nan([numdirs,nsamples]);
    S=zeros([samplesperperiod/2,nsamples]);
    F=nan([numf,nsamples]);
    P=nan([2,nsamples]);
    AF=cell(1,nsamples);
    SOM=zeros(8,nsamples);
    n1=1;
    while n1<=nsamples,
        if ~rem(n1,1e3),disp(n1);
            subplot(121);hist(sum(S(:,1:n1),1),100);
            subplot(122);plot(F(1,1:n1),F(2,1:n1),'.');
            drawnow;
        end
        x=scale*randn(numdirs,1);
        X(:,n1)=x;
        
        [Aud,Som,Outline,af,filt]=diva_synth([X(:,n1);zeros(3,1)],'explicit');
        %x=vt.Scale(idx_art).*x;
        %Outline=vt.Average+vt.Base(:,idx_art)*x;
        
        % computes area function
        %[nill,nill,af,d]=try03_xy2ab(Outline);
        AF{n1}=af;
        SOM(:,n1)=Som;
        placeart=[find(af<=0,1,'first'),find(af<=0,1,'last')];
        if ~isempty(placeart),
            P(:,n1)=placeart;
        else
            P(:,n1)=nan;
        end
        %[filt,f]=try03_a2h(max(min_af,af),d,samplesperperiod,fs);
        filt(1)=0;
        filt=abs(filt(1:end/2));
        S(:,n1)=filt;
        %[filt,f]=try03_a2h0(max(min_af,af),d,samplesperperiod,fs);
        %filt(1)=0;
        %filt=abs(filt(1:end/2));
        
        bfilt=zeros(size(filt)); 
        % JI: next line finds peaks in the filter (formants)
        bfilt(2:end-1)=(filt(2:end-1)>filt(1:end-2)&filt(2:end-1)>filt(3:end));
        idx=find(bfilt);
        if sum(filt)>0&&numel(idx)>=numf,
            f=idx(1:numf)*fs/samplesperperiod;
            F(:,n1)=f;
        else
            F(:,n1)=nan;
        end
        if isempty(placeart),%||(all(af(P(1,n1):P(2,n1))<=0) && (P(2,n1)-P(1,n1))<20 && P(1,n1)>80),
%             subplot(211);plot(Outline);axis equal;
%             subplot(212);plot((0:numel(filt)-1)/numel(filt)*fs/2,filt,'.-');if ~any(isnan(F(:,n1))), set(gca,'xtick',F(:,n1)); end; grid on; 
%             pause;
            n1=n1+1;
        end
    end
    
    k=zeros(1,nsamples);for n1=1:nsamples,k(n1)=numel(AF{n1});end;maxk=max(k);af=zeros(maxk,nsamples);for n1=1:nsamples,af(1:numel(AF{n1}),n1)=AF{n1};af(numel(AF{n1})+1:maxk,n1)=AF{n1}(end);end
    AF=af;
%     save create_mapsdiva01.mat S P F X SOM AF idx_art nsamples scale numf samplesperperiod;
end

%% gm distribution done here - also all of the saving of data
if any(STEPS==2),
    load create_mapsdiva01.mat idx_art nsamples scale numf samplesperperiod S P F X AF idx_art nsamples scale numf samplesperperiod;

    Art=X';
    Fmt=F(1:3,:)';
    idx=find(all(Art>-3&Art<3,2)&~any(isnan(Fmt),2));
    Art=Art(idx,:);
    Fmt=Fmt(idx,:);

    xy=Art;
    sxy=std(xy);
    xy=xy(:,sxy>0);
    ncl=32;
    maxiter=1000;
    options=struct('Display','iter','MaxIter',maxiter);
    gmfit_art=gmdistribution.fit(xy,ncl,'Start','randSample','Replicates',1,'SharedCov',true,'Regularize',1e-5,'Options',options);
    
    save create_mapsdiva01step2.mat gmfit_art 

    x=gmfit_art.posterior(xy); % JI: get posterior probabilities of random inputs
    x=reshape(bsxfun(@times, x, permute([Art,ones(size(Art,1),1)],[1,3,2])),size(x,1),[]);
    B_fmt=pinv(x'*x)*x'*Fmt;
    r_fmt=corrcoef([Fmt,x*B_fmt]);
    
    save create_mapsdiva01step2.mat gmfit_art B_fmt
    
    Som=SOM(1:6,idx)';
    B_som=pinv(x'*x)*x'*Som;
    r_som=corrcoef([Som,x*B_som]);
    
    save create_mapsdiva01step2.mat B_som -append

    if 0,
        iSigma=pinv(gmfit_art.Sigma);MU=gmfit_art.mu;P=gmfit_art.PComponents.';
        Fmt_fit=zeros(size(Fmt));
        for n1=1:size(Art,1),
            dx=bsxfun(@minus,Art(n1,:),MU);
            p=-sum((dx*iSigma).*dx,2)/2;
            p=P.*exp(p-max(p));
            p=p/sum(p);
            px=p*[Art(n1,:),1];
            Fmt_fit(n1,:)=px(:)'*B;
        end
    end
    
    iSigma=pinv(gmfit_art.Sigma);MU=gmfit_art.mu;P=gmfit_art.PComponents.';
    fmfit=struct('mu',MU,'iSigma',iSigma,'p',P,'beta_fmt',B_fmt','beta_som',B_som');
    save diva_synth.mat fmfit -append
end
