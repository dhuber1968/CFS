function FitExp2Data
    
      % assign RT data

    RTs=[2708.622 2317.149 2560.587];
    RTErrors=[149.53232 98.29219 128.4517]; % standard error of the mean difference (Cousinea Morey method)
    RTSD=RTErrors.*sqrt(28); % convert standard error the mean difference to standard deviation of the mean difference

    % read in raw count data and convert to probabilities
    
    D=csvread('RawData.csv');
   
    % convert frequency data to probabilities

    FreqData=zeros(15,4);
    Probs=zeros(28,15,4);
    for i=1:28;  % step through 28 subjects
        for report=1:15
            for cond=1:4
                RawData(i,report,cond)=D(D(:,1)==i & D(:,2)==report & D(:,4)==cond , 3);
            end
        end
        for cond=1:4
            sum(RawData(i,:,cond))
            if cond==1
                Probs(i,:,cond)=RawData(i,:,cond)./sum(RawData(i,:,cond));
            else
                Probs(i,:,cond)=(RawData(i,:,cond)./sum(RawData(i,:,cond)))-Probs(i,:,1);        % differences from No Gestalt condition
            end
        end;
    end
    
    Errors=squeeze(std(Probs,0,1)./sqrt(28));
    FreqData=squeeze(sum(RawData,1));
    
    N=sum(FreqData,1);    
    
    RespOrder=[4 2 12 14 5 6 8 7 3 1 11 13 9 10 15]; % singles, sides, all-but-one, diags, all
    
    params=ones(1,7);
    params=[-0.7128   -2.3084   -1.7053   -2.0994   -3.9151   -1.3087]; % initial parameter values
    Plot=0;                                             % suppress plotting
   [params,err]=fminsearch(@CalcError,params);          % initial fitting
   [params,err]=fminsearch(@CalcError,params);          % second fitting to ensure best fit
   err                                                  % output G-squared error
   params                                               % output best-fitting paramter values
   [(1 - normcdf(0, params(1), 1)) exp(params(2:6)) ]   % convert best-fitting parameters from gammas to probabilities
    Plot=1;                                             % turn on plotting
    CalcError(params);                                  % call error function to produce plotting
    
    function err=CalcError(params)

        %%%%% free parameters %%%%%%
      
        Z_spread =  params(1);    % breakthrough spread to other biters
        SPREAD = 1 - normcdf(0, Z_spread, 1); % convert to probability
        
        Gam_UR =  params(2);        % upper right set to outward pointing
        Gam_UL =  params(2);        % upper left set to outward pointing
        Gam_BR =  params(2);        % bottom right set to outward pointing
        Gam_BL =  params(2);        % bottom left set to outward pointing
        Gam_EYE =  params(3);       % all four biters from eye movements
        Gam_SAL =  params(4);       % salience for inward pointing biter
        Gam_IC  =  params(5);       % both ends of illusory contour (or collinearity)
        Gam_X  =  params(6);        % all biters from illusory X or from eye movements
    
        
        for cond=1:4

            % default values, used for both hate condition
            ALL = exp(Gam_EYE);  % baseline tendency to see all 4 at once
            IRC = 0;  % No illusory left contour (rectangle or collinearity)
            ILC = 0;  % No iIllusory right contour (rectangle or collinearity)
            UR = exp(Gam_UR); % upper right biter
            UL = exp(Gam_UL); % upper left biter
            BR = exp(Gam_BR); % bottom right biter
            BL = exp(Gam_BL); % bottom left biter

            if  cond==2 % both talk condition
            %    ILC = exp(Gam_IC); % include for model with collinearity
            %    IRC = exp(Gam_IC); % include for model with collinearity

                ALL= exp(Gam_X);   % include for model with illusory X

        % include the following for salience model
                UR = exp(Gam_SAL);  % increase salience of inward pointing biters
                UL = exp(Gam_SAL);
                BR = exp(Gam_SAL);
                BL = exp(Gam_SAL);
            elseif cond==3 % left talk condition
                ILC = exp(Gam_IC);      % include for collinearity or illusory rectangle

        % include the following for salience model
                UL = exp(Gam_SAL); % increase salience of inward pointing biters
                BR = exp(Gam_SAL);
            elseif cond==4 % right talk condition
                IRC = exp(Gam_IC);       % include for collinearity or illusory rectangle

        % include the following for salience model
                UR = exp(Gam_SAL);  % increase salience of inward pointing biters
                BL = exp(Gam_SAL);
            end
            [TotResp RT]=CalcResp(ALL,IRC,ILC,UR,UL,BR,BL,SPREAD);
            
            if cond<4 % condition 3 and 4 make the same RT prediction
                RTPred(cond)=RT;
            end
            
            FreqPred(:,cond)=N(cond).*TotResp;
                     
            ProbPred(:,cond)=FreqPred(:,cond)./N(cond);
            ProbData(:,cond)=FreqData(:,cond)./N(cond);
        end
                
        if Plot==1
            CondOrder=[1 3 4 2];
            figure(1);
            for c=1:4
                cond=CondOrder(c);
                subplot(4,1,c);
                hold off
                if cond==1
                    plot(FreqPred(RespOrder,cond)./N(cond),'xr','LineWidth',2,'MarkerSize',10);
                else
                    plot(FreqPred(RespOrder,cond)./N(cond)-FreqPred(RespOrder,1)./N(1),'xr','LineWidth',2,'MarkerSize',10);
                end
                hold on
                
                set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);
                if cond==1
                    plot(FreqData(RespOrder,cond)./N(cond),'ob','LineWidth', 2);
                    errorbar([1:15],FreqData(RespOrder,cond)./N(cond),Errors(:,cond),Errors(:,cond),'LineStyle','none','LineWidth', 2,'Color','b');
                else
                    plot(FreqData(RespOrder,cond)./N(cond)-FreqData(RespOrder,1)./N(1),'ob','LineWidth', 2);
                    errorbar([1:15],FreqData(RespOrder,cond)./N(cond)-FreqData(RespOrder,1)./N(1),Errors(:,cond),Errors(:,cond),'LineStyle','none','LineWidth', 2,'Color','b');
                end
                
                legend('Model','Data','Location','NorthWest');
                if cond==1
                    axis([.5 15.5 0 .4]);
                    line([4.5 4.5],[0 .4]);
                    line([8.5 8.5],[0 .4]);
                    line([12.5 12.5],[0 .4]);
                else
                    axis([.5 15.5 -.1 .1]);
                    line([.5 15.5],[0 0]);
                    line([4.5 4.5],[-.1 .1]);
                    line([8.5 8.5],[-.1 .1]);
                    line([12.5 12.5],[-.1 .1]);
                end
                if cond==1
                    title('No Gestalt Condition');
                elseif cond==2
                    title('Illusory Cross Condition');
                elseif cond==3
                    title('Illusory Left Rectangle Condition');
                elseif cond==4
                    title('Illusory Right Rectangle Condition');
                end
                set(gca,'XTick',[1:15]);
                 ax = gca;
                 ax.XTick = [1:15];
                 ax.XTickLabel = '';
                myLabels = {'+ -','- +','- -','- -','+ +','- +','- -','+ -','- +','+ -','+ +','+ +','+ -','- +','+ +';
                            '- -','- -','- +','+ -','- -','- +','+ +','+ -','+ +','+ +','+ -','- +','- +','+ -','+ +';
                            '','One Biter Reported','','','','One Side Reported','','','','All But One Reported','','','','Gestalt Reported',''};
                        
                for i=1:15
                    text(i, ax.YLim(1), sprintf('%s\n%s\n%s', myLabels{:,i}), ...
                            'horizontalalignment', 'center', 'verticalalignment', 'top','FontName','Arial','FontSize',12,'FontWeight','Bold');  
                end

                if cond==1
                    ylabel('Report Probability');
                else
                    ylabel('Change from No Gestalt Condition');
                end
                
                filename=sprintf('Cond%d.bmp',cond);
                axes('pos',[.45 (.85-(c-1)*.22) .08 .08])
                imshow(filename);
                
                
            end
            figure(2);
            hold off
            plot(RTPred([1 3 2]),'xr','LineWidth',2,'MarkerSize',10);
            hold on
            
            set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);
            plot(RTs([1 3 2]),'ob','LineWidth', 2);
            errorbar([1:3],RTs([1 3 2]),RTErrors([1 3 2]),RTErrors([1 3 2]),'LineStyle','none','LineWidth', 2,'Color','b');
            axis([.5 3.5 2100 2900]);
            legend('Model','Data','Location','SouthWest');
             ax = gca;
             ax.XTick = [1:3];
             ax.XTickLabel = ({'No Gestalt' 'Illusory Rectangle' 'Illusory Cross'});
             ax.YTick = [2200:200:3000];
             ylabel('Median Breakthrough Latency (ms)');
             xlabel('Condition');
            
             % if plotting, also output percent variance accounted for
            tot_mean=mean(ProbData(:));
            SStot=sum((ProbData(:)-mean(ProbData(:))).^2);
            SSres=sum((ProbData(:)-ProbPred(:)).^2);
            ProbExplained=1-(SSres/SStot)
            
            SStot=sum((RTs(:)-mean(RTs(:))).^2);
            SSres=sum((RTs(:)-RTPred(:)).^2);
            RTExplained=1-(SSres/SStot)
        end
        % calculate G-squared error based on likelihood ratio test
        err=0;
        err=-2*sum(log(mnpdf(FreqData',ProbPred')./mnpdf(FreqData',ProbData')));
        err=err-2*(28*70)*sum(log(normpdf(RTPred(:),RTs(:),RTSD(:))./normpdf(RTs(:),RTs(:),RTSD(:))));
        % each of the 3 RT datapoints is from 28 subjects and 70 trials so
        % the likelihood is multiplied by (28*70) to be comparable to the
        % likelihood of the raw response data in terms of trials
    end

    function [TotResp RT]=CalcResp(ALL,IRC,ILC,UR,UL,BR,BL,SPREAD)
      
        % calculate chances of spread occuring (N) or not (n) to each of
        % the other 2 or 3 biters
        nnn = (1-SPREAD).^3;           % none of other three
        Nnn = SPREAD*((1-SPREAD).^2);   
        NNn = (SPREAD.^2)*(1-SPREAD);    
        NNN = SPREAD.^3;               % spread to all of the other three
        nn = (1-SPREAD).^2;   
        Nn = SPREAD*(1-SPREAD);  
        NN = SPREAD.^2;      
        
        NONE=1;  % perceptual salience of Mondrian
        
        SamplingSpace=NONE+ALL+IRC+ILC+UR+UL+BR+BL; % denominator for softmax
        
        NONE=NONE/SamplingSpace;
        ALL=ALL/SamplingSpace;
        IRC=IRC/SamplingSpace;
        ILC=ILC/SamplingSpace;
        UR=UR/SamplingSpace;
        UL=UL/SamplingSpace;
        BR=BR/SamplingSpace;
        BL=BL/SamplingSpace;

        % perception probabilities for each cycle
        
        resp(2) = UR*nnn;                                          % upperR only
        resp(4) = UL*nnn;                                          % upperL only
        resp(12)= BR*nnn;                                          % bottomR only
        resp(14)= BL*nnn;                                          % botomL only
        resp(5) = (UL+UR)*Nnn;                                 % both top biters
        resp(6) = (UR+BR)*Nnn;                                 % both right biters
        resp(7) = (UL+BL)*Nnn;                                 % both left biters
        resp(8) = (BL+BR)*Nnn;                                 % both bottom biters
        resp(3) = IRC*Nn + (UR+BL+BR)*NNn;                        % all but upper left
        resp(1) = ILC*Nn +  (UL+BR+BL)*NNn;                        % all but upper right 
        resp(11)= IRC*Nn + (UL+UR+BL)*NNn;                        % all but bottom right
        resp(13)= ILC*Nn + (UL+UR+BR)*NNn;                         % all but bottom left
        resp(9) = ILC*nn + (UL+BR)*Nnn;                           % diagL
        resp(10) = IRC*nn + (UR+BL)*Nnn;                          % diagR
        resp(15)= (ILC+IRC)*NN + ALL + NNN*(UR+UL+BR+BL);         % all biters

        if abs((1-sum(resp))-NONE)>.00000001  % check to see if sum of 1-15 = none
            disp('calculation error');
        end      
        
        TotResp=resp./sum(resp); % because cycling is done until breakthrough, the final
        % response distrbitution is the normalized response distribution,
        % excluding the NONE
   
        RT = 1./sum(resp); % cycling til breakthrough is a geometric distribution and
                    %  1/B is the mean, where B is the probability of breakthrough on each cycle

        RT = RT *1000; % convert to milliseconds
    end

end
