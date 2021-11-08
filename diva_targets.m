function varargout=diva_targets(option,varargin)
global DIVA_x;

switch(lower(option))
    case 'list' % target files in diva_targets folder
        filename=fullfile(fileparts(which(DIVA_x.model)),[DIVA_x.model,'.csv']);
        [production_ids,production_labels]=diva_targets_readcsvfile(filename);
        varargout{1}=production_labels;
        varargout{2}=production_ids;
    case 'delete'
        filename=fullfile(fileparts(which(DIVA_x.model)),[DIVA_x.model,'.csv']);
        [production_ids,production_labels]=diva_targets_readcsvfile(filename);
        production=varargin{1};
        idx=strmatch(production,production_labels,'exact');
        if ~isempty(idx),
            idx=idx(1);
            if isnan(production_ids(idx)), 
                filename_mat1=fullfile(fileparts(which(DIVA_x.model)),[DIVA_x.model,'_targets'],production_labels{idx});
                filename_mat2=fullfile(fileparts(which(DIVA_x.model)),[DIVA_x.model,'_targets'],['bak_',production_labels{idx}]);
                if ~isdir(fileparts(filename_mat1)), [nill,nill]=mkdir(fileparts(filename_mat1)); end
            else 
                filename_mat1=fullfile(fileparts(which(DIVA_x.model)),[DIVA_x.model,'_',num2str(production_ids(idx),'%06d')]);
                filename_mat2=fullfile(fileparts(which(DIVA_x.model)),['bak_',DIVA_x.model,'_',num2str(production_ids(idx),'%06d')]);
            end
            production_ids=production_ids([1:idx-1,idx+1:numel(production_ids)]);
            production_labels=production_labels([1:idx-1,idx+1:numel(production_labels)]);
            diva_targets_writecsvfile(filename,production_ids,production_labels);
            if isunix
                [nill,ok]=system(['mv ',filename_mat1,'.mat ',filename_mat2,'.mat']); if nill,disp(ok); end
                [nill,ok]=system(['mv ',filename_mat1,'.txt ',filename_mat2,'.txt']); if nill,disp(ok); end
                fprintf('Target %s deleted (backup stored as %s)\n',filename_mat1,filename_mat2);

            else
                [nill,ok]=system(['move ',filename_mat1,'.mat ',filename_mat2,'.mat']); if nill,disp(ok); end
                [nill,ok]=system(['move ',filename_mat1,'.txt ',filename_mat2,'.txt']); if nill,disp(ok); end
                fprintf('Target %s deleted (backup stored as %s)\n',filename_mat1,filename_mat2);
            end
        else
            disp(['warning: no match for production ',production,' in ',filename]);
        end
    case 'new'
        filetype=varargin{1};
        switch(lower(filetype))
            case 'txt'
                varargout{1}=diva_targets_initstruct;
            case 'mat'
                production_info=diva_targets_initstruct;
                varargout{1}=diva_targets('timeseries',production_info,'header');
        end
    case 'random'
        filetype=varargin{1};
        switch(lower(filetype))
            case 'txt'
                varargout{1}=diva_targets_initstruct(3);
            case 'mat'
                production_info=diva_targets_initstruct(3);
                varargout{1}=diva_targets('timeseries',production_info,'header');
        end
    case 'load'
        filename=fullfile(fileparts(which(DIVA_x.model)),[DIVA_x.model,'.csv']);
        [production_ids,production_labels]=diva_targets_readcsvfile(filename);
        filetype=varargin{1};
        production=varargin{2};
        idx=strmatch(production,production_labels,'exact');
        if ~isempty(idx)
            if numel(idx)>1, disp(['warning: multiple entries matching ',production,' in ',filename]); idx=idx(1); end
            if isnan(production_ids(idx)),
                filename=fullfile(fileparts(which(DIVA_x.model)),[DIVA_x.model,'_targets'],production_labels{idx});
                if ~isdir(fileparts(filename)), [nill,nill]=mkdir(fileparts(filename)); end
           else
                filename=fullfile(fileparts(which(DIVA_x.model)),[DIVA_x.model,'_',num2str(production_ids(idx),'%06d')]);
            end
            switch(lower(filetype))
                case 'txt'
                    filename=[filename,'.txt'];
                    if ~isempty(dir(filename))
                        varargout{1}=diva_targets_txt2struct(filename);
                    else
                        varargout{1}=[];
                    end
                case 'mat'
                    filename=[filename,'.mat'];
                    if ~isempty(dir(filename))
                        load(filename,'timeseries');
                        varargout{1}=timeseries;
                    else
                        production_info=diva_targets_txt2struct(filename);
                        timeseries=diva_targets('timeseries',production_info,'header');
                        varargout{1}=timeseries;
                    end
            end
        else
            disp(['warning: no entry matching ',production,' in ',filename]); 
        end
    case 'reset'
        production_info=varargin{1};
        varargout={diva_targets('timeseries',production_info,'header')};
    case 'save'
        filename=fullfile(fileparts(which(DIVA_x.model)),[DIVA_x.model,'.csv']);
        [production_ids,production_labels]=diva_targets_readcsvfile(filename);
        filetype=varargin{1};
        production=varargin{2};
        production_info=varargin{3};
        if nargin<5, overwrite=-1; else overwrite=varargin{4}; end
        idx=strmatch(production,production_labels,'exact');
        if numel(idx)>1, disp(['warning: multiple entries matching ',production,' in ',filename]); idx=idx(1); end
        if isempty(idx)
            %if isempty(production_ids), production_ids=1; else production_ids(end+1)=max(production_ids)+1; end
            if isempty(production_ids), production_ids=nan; else production_ids(end+1)=nan; end % note: default new targets in diva_targets
            production_labels(end+1)={production};
            diva_targets_writecsvfile(filename,production_ids,production_labels);
            idx=numel(production_ids);
        end
        if isnan(production_ids(idx)),
            filename=fullfile(fileparts(which(DIVA_x.model)),[DIVA_x.model,'_targets'],production_labels{idx});
            if ~isdir(fileparts(filename)), [nill,nill]=mkdir(fileparts(filename)); end
        else
            filename=fullfile(fileparts(which(DIVA_x.model)),[DIVA_x.model,'_',num2str(production_ids(idx),'%06d')]);
        end
        
        switch(lower(filetype))
            case 'txt'
                production_info.name=production;
                filename=[filename,'.txt'];
                diva_targets_struct2txt(production_info,filename);
            case 'mat'
                production_info.name=production;
                filename=[filename,'.mat'];
                timeseries=diva_targets('timeseries',production_info,'header');
                if ~isempty(dir(filename)) && overwrite~=1
                    old=load(filename,'-mat');
                    if ~overwrite
                        if size(old.timeseries.Art,1)==size(timeseries.Art,1), 
                            timeseries.Art=old.timeseries.Art;
                        else
                            disp('warning: incorrect length of learned Feedforward timeseries. overwriting');
                        end
                    elseif size(old.timeseries.Art,1)==size(timeseries.Art,1), 
                        ok=questdlg('Overwrite existing Articulatory sequence?','','Yes', 'No', 'No');
                        if strcmp(ok,'No')
                            timeseries.Art=old.timeseries.Art;
                        end
                    end
                end
                save(filename,'timeseries');
            case 'art'
                filename=[filename, '.mat'];
                if ~isempty(dir(filename))
                    load(filename,'timeseries');
                    timeseries.Art=production_info;
                    save(filename,'timeseries');
                elseif numel(varargin)>=4
                    timeseries=diva_targets('timeseries',varargin{4},'header');
                    timeseries.Art=production_info;
                    save(filename,'timeseries');
                else
                    disp(['warning: file ',filename,' not found']);
                end
        end
    case 'timeseries'
        production_info=varargin{1};
        if nargin<3, doheader=0; else doheader=strcmpi(varargin{2},'header'); end
        params=diva_vocaltract;
        [Aud_min,Aud_max,Time]=diva_targets_timeseriesconvert(production_info,params.Output(1));
        [Som_min,Som_max,Time]=diva_targets_timeseriesconvert(production_info,params.Output(2));
        Art=repmat(mean(params.Input.Default,2)',[size(Aud_min,1),1]);
        N_samplesperheader=0;
        if doheader&&N_samplesperheader>0
            k0=1-(1-linspace(0,1,N_samplesperheader)').^2;
            Aud_min=cat(1,k0*Aud_min(1,:)+(1-k0)*params.Output(1).Default(:,1)', Aud_min, flipud(k0)*Aud_min(end,:)+(1-flipud(k0))*params.Output(1).Default(:,1)');
            Aud_max=cat(1,k0*Aud_max(1,:)+(1-k0)*params.Output(1).Default(:,2)', Aud_max, flipud(k0)*Aud_max(end,:)+(1-flipud(k0))*params.Output(1).Default(:,2)');
            Som_min=cat(1,k0*Som_min(1,:)+(1-k0)*params.Output(2).Default(:,1)', Som_min, flipud(k0)*Som_min(end,:)+(1-flipud(k0))*params.Output(2).Default(:,1)');
            Som_max=cat(1,k0*Som_max(1,:)+(1-k0)*params.Output(2).Default(:,2)', Som_max, flipud(k0)*Som_max(end,:)+(1-flipud(k0))*params.Output(2).Default(:,2)');
            Art=cat(1,Art(1+zeros(N_samplesperheader,1),:),Art,Art(end+zeros(N_samplesperheader,1),:));
            dt=min(diff(Time));
            Time=cat(1,dt*(0:N_samplesperheader-1)',dt*N_samplesperheader+Time,dt*N_samplesperheader+Time(end)+dt*(1:N_samplesperheader)');
        end
        varargout{1}=struct('Aud_min',Aud_min,'Aud_max',Aud_max,'Som_min',Som_min,'Som_max',Som_max,'Art',Art,'time',Time);
    case 'combine'
        production_info1=varargin{1};
        production_info2=varargin{2};
        if numel(varargin)<3, duration_transition=40; else duration_transition=varargin{3}; end
        if numel(varargin)<4, duration2=[]; else duration2=varargin{4}; end
        if numel(varargin)<5, timeseries1=[]; else timeseries1=varargin{5}; end
        if numel(varargin)<6, timeseries2=[]; else timeseries2=varargin{6}; end
        %production_info <- production_info1 & production_info2;
        production_info.name='';
        if ~isempty(duration2), production_info2.length=duration2; end
        production_info.length=production_info1.length+duration_transition+production_info2.length;
        production_info.wrapper=0;
        production_info.interpolation='linear';
        if ~isequal(production_info1.interpolation,'linear'), fprintf('warning: interpolation methods changed from %s to linear\n',production_info1.interpolation); end
        if ~isequal(production_info2.interpolation,'linear'), fprintf('warning: interpolation methods changed from %s to linear\n',production_info2.interpolation); end
        params=diva_vocaltract;
        for n0=1:numel(params.Output),
            for n1=1:numel(params.Output(n0).Plots_dim)
                if numel(params.Output(n0).Plots_dim{n1})==1
                    idx=params.Output(n0).Plots_dim{n1};
                    x0=production_info1.([params.Output(n0).Plots_label{n1},'_control']);
                    x1=production_info1.([params.Output(n0).Plots_label{n1},'_min']);
                    x2=production_info1.([params.Output(n0).Plots_label{n1},'_max']);
                    if isequal(x1,'X'), x1=params.Output(n0).DefaultSound(idx,1); end
                    if isequal(x2,'X'), x2=params.Output(n0).DefaultSound(idx,2); end
                    if any(isnan(x1)), x1(isnan(x1))=params.Output(n0).DefaultSound(idx,1); end
                    if any(isnan(x2)), x2(isnan(x2))=params.Output(n0).DefaultSound(idx,2); end
                    y0=production_info2.([params.Output(n0).Plots_label{n1},'_control']);
                    y1=production_info2.([params.Output(n0).Plots_label{n1},'_min']);
                    y2=production_info2.([params.Output(n0).Plots_label{n1},'_max']);
                    if isequal(y1,'X'), y1=params.Output(n0).DefaultSound(idx,1); end
                    if isequal(y2,'X'), y2=params.Output(n0).DefaultSound(idx,2); end
                    if any(isnan(y1)), y1(isnan(y1))=params.Output(n0).DefaultSound(idx,1); end
                    if any(isnan(y2)), y2(isnan(y2))=params.Output(n0).DefaultSound(idx,2); end
                    if numel(x0)==1&&numel(y0)==1&&x1==y1&&x2==y2,
                    else
                        if numel(x0)==1, x0=[0 production_info1.length]; x1=[x1 x1]; x2=[x2 x2]; end
                        if numel(y0)==1, y0=[0 production_info2.length]; y1=[y1 y1]; y2=[y2 y2]; end
                        x0=[x0 production_info1.length+duration_transition+y0];
                        x1=[x1 y1];
                        x2=[x2 y2];
                    end
                    production_info.([params.Output(n0).Plots_label{n1},'_control'])=x0;
                    production_info.([params.Output(n0).Plots_label{n1},'_min'])=x1;
                    production_info.([params.Output(n0).Plots_label{n1},'_max'])=x2;
                end
            end
        end
        %timeseries <- timeseries1 & timeseries2;
        timeseries=diva_targets('timeseries',production_info,'header');
        if ~isempty(timeseries1),
            t0=timeseries.time<=production_info1.length;
            timeseries.Art(t0,:)=timeseries1.Art(min(size(timeseries1.Art,1),1:nnz(t0)),:);
        end
        if ~isempty(timeseries2),
            t0=timeseries.time>=production_info1.length+duration_transition;
            timeseries.Art(t0,:)=timeseries2.Art(min(size(timeseries2.Art,1),1:nnz(t0)),:);
        end
        t0=find(timeseries.time>=production_info1.length & timeseries.time<=production_info1.length+duration_transition);
        if t0(1)>1, t0=[t0(1)-1;t0(:)]; end
        if t0(end)<size(timeseries.Art,1), t0=[t0(:);t0(end)+1]; end
        for n1=1:size(timeseries.Art,2), timeseries.Art(t0,n1)=linspace(timeseries.Art(t0(1),n1),timeseries.Art(t0(end),n1),numel(t0))'; end
        varargout={production_info, timeseries};
    otherwise
        varargout{1}=feval([mfilename,'_',option],varargin{:});
end
end


function production_info=diva_targets_initstruct(nrandom)
if nargin<1, nrandom=0; end
production_info.name='';
production_info.length=500;
production_info.wrapper=0;
production_info.interpolation='spline';
params=diva_vocaltract;
for n0=1:numel(params.Output),
    for n1=1:numel(params.Output(n0).Plots_dim)
        if numel(params.Output(n0).Plots_dim{n1})==1
            idx=params.Output(n0).Plots_dim{n1};
            if isfield(params.Output(n0),'DefaultSound'), Default=params.Output(n0).DefaultSound; else Default=params.Output(n0).Range; end
            if n0==1&&nrandom>0
                production_info.([params.Output(n0).Plots_label{n1},'_control'])=linspace(0,production_info.length,nrandom);
                x=sort(rand(2,nrandom));
                production_info.([params.Output(n0).Plots_label{n1},'_min'])=Default(idx,1)*(1-x(1,:))+Default(idx,2)*x(1,:);
                production_info.([params.Output(n0).Plots_label{n1},'_max'])=Default(idx,1)*(1-x(2,:))+Default(idx,2)*x(2,:);
            else
                production_info.([params.Output(n0).Plots_label{n1},'_control'])=0;
                production_info.([params.Output(n0).Plots_label{n1},'_min'])=Default(idx,1);
                production_info.([params.Output(n0).Plots_label{n1},'_max'])=Default(idx,2);
            end
        end
    end
end
end

function [y_min,y_max,Time]=diva_targets_timeseriesconvert(production_info,params_info)
DT=5;
Nt=1+ceil((production_info.length+2*production_info.wrapper)/DT);
Time=(0:Nt-1)'*DT;
y_min=zeros([Nt,params_info.Dimensions]);
y_max=zeros([Nt,params_info.Dimensions]);
for n1=1:numel(params_info.Plots_dim)
    if numel(params_info.Plots_dim{n1})==1
        idx=params_info.Plots_dim{n1};
        x0=production_info.([params_info.Plots_label{n1},'_control'])+production_info.wrapper;
        x1=production_info.([params_info.Plots_label{n1},'_min']);
        x2=production_info.([params_info.Plots_label{n1},'_max']);
        if isequal(x1,'X'), x1=params_info.DefaultSound(idx,1); end
        if isequal(x2,'X'), x2=params_info.DefaultSound(idx,2); end
        if any(isnan(x1)), x1(isnan(x1))=params_info.DefaultSound(idx,1); end
        if any(isnan(x2)), x2(isnan(x2))=params_info.DefaultSound(idx,2); end
        y_min(:,idx)=diva_targets_interpolate(x0,x1,Time,production_info.interpolation);
        y_max(:,idx)=diva_targets_interpolate(x0,x2,Time,production_info.interpolation);
    end
end
if production_info.wrapper>0, 
    mask=Time<production_info.wrapper | Time>production_info.wrapper+production_info.length;
    y0_min=params_info.DefaultSilence(:,1);
    y0_max=params_info.DefaultSilence(:,2);
    y_min(mask,:)=repmat(y0_min',nnz(mask),1);
    y_max(mask,:)=repmat(y0_max',nnz(mask),1);
end
temp=max(y_min,y_max);
y_min=min(y_min,y_max);
y_max=temp;
end

function y=diva_targets_interpolate(x0,y0,x,interpolation)
x0=x0(:);
y0=y0(:);
x=x(:);
if numel(x0)~=numel(y0), y=nan(size(x)); return; end
[x0,idx]=sort(x0);
y0=y0(idx);
idx=find(x0(2:end)==x0(1:end-1));
x0(idx)=x0(idx)+eps;
if min(x0)>0,x0=[0;x0];y0=[y0(1);y0]; end
if max(x0)<max(x),x0=[x0;max(x)];y0=[y0;y0(end)]; end
if numel(x0)<=2, interpolation='linear'; end
if numel(x0)<=1, interpolation='nearest'; end
y=interp1(x0,y0,x,interpolation);
end

function [production_id,production_label]=diva_targets_readcsvfile(filename)
filepath=regexprep(filename,'\.csv$','_targets');
production_id=[];
production_label={};
if isdir(filepath) % targets in diva_targets folder
    fname=dir(fullfile(filepath,'*.txt'));
    [nill,tname,nill]=cellfun(@fileparts,{fname.name}','uni',0);
    tname=tname(cellfun('length',regexp(tname,'^bak_'))==0);
    production_label=[production_label;tname(:)];
    production_id=[production_id; nan(numel(tname),1)];
end
if ~isempty(dir(filename)) % targets listed in diva.csv file
    [tproduction_id,tproduction_label]=textread(filename,'%n%s','delimiter',',','headerlines',1);
    production_label=[production_label;tproduction_label(:)];
    production_id=[production_id; tproduction_id];
elseif isempty(production_id)
    disp(['warning: targets file ',filename,' or directory ',filepath,' does not exist or is empty: initializing']);
end
end

function diva_targets_writecsvfile(filename,production_id,production_label)
if any(~isnan(production_id))
    fh=fopen(filename,'wt');
    fprintf(fh,'ID,Label\n');
    for n1=1:numel(production_id)
        if ~isnan(n1), fprintf(fh,'%d,%s\n',production_id(n1),production_label{n1}); end
    end
    fclose(fh);
end
end

function out=diva_targets_txt2struct(filename,out)
if nargin<2, out=[]; end
comment=0;
fieldname='arg';
s=textread(filename,'%s');
for n1=1:length(s),
    if comment || isempty(s{n1}),
    elseif strncmp(s{n1},'%{',2), % comment open
        comment=1;
    elseif strncmp(s{n1},'%}',2), % comment close
        comment=0;
    elseif s{n1}(1)=='#', % field name
        fieldname=(s{n1}(2:end));
        out.(fieldname)=[];
    else % field value
        n=str2double(s{n1});
        if ~isnan(n)&&all(ismember(s{n1},'0123456789.+-')), newvalue=n; else newvalue=s{n1}; end; %{s{n1}}; end
        if isfield(out,fieldname) && ~isempty(out.(fieldname)),
            out.(fieldname)=cat(2,out.(fieldname),newvalue);
        else
            out.(fieldname)=newvalue;
        end
    end
end
end

function diva_targets_struct2txt(out,filename)
fh=fopen(filename,'wt');
s=fieldnames(out);
for n1=1:length(s),
    fprintf(fh,'#%s\n',s{n1});
    if ischar(out.(s{n1}))
        fprintf(fh,'%s\n',out.(s{n1}));
    else
        x=out.(s{n1});
        for n2=1:numel(x),
            n=sum(rem(x(n2)*logspace(0,6,7),1)~=0);
            fprintf(fh,['%0.',num2str(n),'f '],x(n2));
        end
        fprintf(fh,'\n');
    end
end
fclose(fh);
end



    

