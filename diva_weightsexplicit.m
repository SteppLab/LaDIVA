function diva_weightsexplicit(block)
% Level-2 M file S-Function.
  setup(block);  
end

%% Initialization   
function setup(block)
  
  % Register number of dialog parameters
  block.NumDialogPrms = 3;
  block.DialogPrmsTunable = {'Nontunable','Nontunable','Nontunable'};

  % Register number of input and output ports
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 1;

  % Setup functional port properties to dynamically inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = -1;
  block.InputPort(1).DirectFeedthrough = true;  
  block.InputPort(2).Dimensions        = -1;
  block.InputPort(2).DirectFeedthrough = true;  
  block.OutputPort(1).Dimensions       = -1;
  
  % Set block sample time to discrete
  block.SampleTimes = [-1 0];
  
  % Register methods
  block.RegBlockMethod('SetInputPortDimensions',  @SetInputDims);
  block.RegBlockMethod('Outputs',                 @Output);  
  
end

 
function SetInputDims(block, port, dm)
    block.InputPort(port).Dimensions = dm;
    if port==2, block.OutputPort(1).Dimensions = dm; end
end

%% Output & Update equations   
function Output(block)

    dy=block.InputPort(1).Data;
    x=block.InputPort(2).Data;
    nout=block.DialogPrm(1).Data;
%     switch(nout)
%         case 1, nout='auditory';
%         case 2, nout='somatosensory';
%     end
%     switch(lower(nout)),
%         case 'auditory', nout=1;
%         case 'somatosensory', nout=2;
%     end
    EPS=block.DialogPrm(2).Data;
    LAMBDA=block.DialogPrm(3).Data;
    N=numel(x);
    M=numel(dy);

    Ix=eye(N);
    Iy=eye(M);
    Q=diva_vocaltract('base');%Ix;%orth(randn(N));
    %Q(13,13) = 1000;
    DX=EPS*Q; % direction of articulatory change
    
    DY=zeros([M,N]); % direction of auditory/somatosensory change
    y=diva_vocaltract(nout,x); % HW an actual motor representation : cannot be negative, so negative should be eleminated
    
%     % set x dimension to zero if x values for tension reach corner cases HRW
%     if x(11,1) == 0 DX(:,11) = zeros(N,1); end
%     if x(12,1) == 0 DX(:,12) = zeros(N,1); end
%     if x(13,1) == 0 DX(:,13) = zeros(N,1); end
    
    for ndim=1:N, % computes jacobian
        xt=x+DX(:,ndim);
       
%         if ndim == 11 || ndim == 12 % remove tuning of CT and TA
%            xt=x;                    % remove tuning of CT and TA
%         end                         % remove tuning of CT and TA
%         
        DY(:,ndim)=diva_vocaltract(nout,xt)-y;
    end
    JJ=DY*DY';
    if sum(sum(isnan(JJ+LAMBDA*EPS^2*Iy)| isinf(JJ+LAMBDA*EPS^2*Iy)))
        disp('jacobian product has NaN or Inf');
    end
    iJ=EPS*Q*DY'*pinv(JJ+LAMBDA*EPS^2*Iy); % computes pseudoinverse
    dx=-iJ*dy;

    
   %   if strcmp(nout, 'Auditory')
% %       %HW : dx and dy are error signals - consider the aboslute value of change when
% %       %considering forward transform.
% %       iJ(11,1) =DY(1,11)*EPS*(1/(JJ(1,1)+LAMBDA*EPS^2));
% %       dx(11,1) = -iJ(11,1)*dy(1,1);
       % disp(['Auditory : DY ',num2str(DY(1,11)),' dy ',num2str(dy(1,1)),' dx ',num2str(dx(11,1)),' y ',num2str(y(1,1)),' x ',num2str(x(11,1))]);
    %  end
%      
%      if strcmp(nout, 'Somatosensory') 
%          %dx(13,1)=1200;
%        disp(['Somatosensory : DY ',num2str(DY(7,13)),' dy ',num2str(dy(7,1)),' dx ',num2str(dx(13,1)),' y ',num2str(y(7,1)),' x ',num2str(x(13,1))]);
%      end
    
    block.OutputPort(1).Data = dx;
%     K=.05;
%     dx0=-K*(Ix-iJ*DY*Q'/EPS)*x;
%     block.OutputPort(1).Data = dx+dx0;

end

