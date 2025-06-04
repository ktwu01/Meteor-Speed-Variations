% load original data
clear all;clc;close all;
load('D:\000learning\00research\00DataAndCode\Mengcheng_wkt_Nhindley\OptionalMMRdata2014_2022.mat');

% Mention! BEFORE YOU RUN THIS CODE
% please turn on ALL or SOME of the if (1) so that you can satisfile your needs!
% try not to just directely run this file instead.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the PCA slope of each day of 9 years
% time use: about 7 mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0)
PCAslope9yrs=zeros(3287,1);
    for i = 1:3287
        id = find(mcest >=24*(i-1) & mcest <24*i);
        if length(id)>1500 % more than 1500 per day will be included
            mch_id=mch(id);      speed_id=speed(id)./1000;
            % limit the data range and let
            % the length of X(:,1) and X(:,2) be the same
            fd= find(mch_id> 70 & mch_id < 110 ...
                & speed_id >8 & speed_id< 80);
            mch_id=mch_id(fd);
            speed_id=speed_id(fd);

            % X is the original data for PCA method
            X=zeros(20000,2);
            X(1:length(fd),1)=speed_id;
            X(1:length(fd),2)=mch_id;
%             n=1;
%             for sp=8:4:80
%                 m=1;
%                 for ht=70:2:110
%                     numdis_id(n,m)=length(find(mch_id>ht-1 & ...
%                         mch_id<=ht+1 & speed_id>sp-2 & speed_id<sp+2));
%                     m=m+1;
%                 end
%                 n=n+1;
%             end
%             
%             % numdis_id=numdis_id'
%             numdis_id=log(numdis_id/max(max(numdis_id)));
%             H1=scatter(X(:,1),X(:,2),'x','LineWidth',2,varargin{:},'visible','off');

% Calcu Main direction
            Xmean=mean(X);
            X0=X-Xmean;
            covMat=(X0.')*X0./size(X,1);
            [V,~]=eigs(covMat,1);
            V=V.*sign(V(1));
            Vc=[V(2);-V(1)];
            Lp=max(X0*V);
            Ln=min(X0*V);
            Ll=max(X0*Vc);
            Lr=min(X0*Vc);
% Draw Histogram
            H2=plot(Xmean(1)+[Lp,Ln].*V(1).*1.05,Xmean(2)+[Lp,Ln].*V(2).*1.05,'visible','off');

%% calcu the slope
            PCAslope9yrs(i)=(max(H2.YData)-min(H2.YData))/(max(H2.XData)-min(H2.XData));
%             % this is the plot of main direction
%             % hold on; plot(H2.XData, H2.YData)

%             % hold on;
%             % H3=fill(Xmean(1)+[Lp,Lp,Ln,Ln].*V(1).*1.05+[Ll,Lr,Lr,Ll].*Vc(1).*1.05,...
%             %      Xmean(2)+[Lp,Lp,Ln,Ln].*V(2).*1.05+[Ll,Lr,Lr,Ll].*Vc(2).*1.05,... 
%             %     H1.CData,'FaceAlpha',.15,'EdgeColor','none');
%             
%             histHdl=histogram(X0*Vc,'BinEdges',linspace(Lr,Ll,13),'Visible','off');
%             maxValue=max(histHdl.Values);
% 
%             for i=1:histHdl.NumBins
%                 H4(i)=fill(Xmean(1)+[histHdl.BinEdges(i).*Vc(1)+Lp.*V(1).*1.05,...
%                     histHdl.BinEdges(i).*Vc(1)+Lp.*V(1).*1.05+histHdl.Values(i).*Lp.*V(1).*0.5./maxValue,...
%                     histHdl.BinEdges(i+1).*Vc(1)+Lp.*V(1).*1.05+histHdl.Values(i).*Lp.*V(1).*0.5./maxValue,...
%                     histHdl.BinEdges(i+1).*Vc(1)+Lp.*V(1).*1.05]./1.5,...
%                     Xmean(2)+[histHdl.BinEdges(i).*Vc(2)+Lp.*V(2).*1.05,...
%                     histHdl.BinEdges(i).*Vc(2)+Lp.*V(2).*1.05+histHdl.Values(i).*Lp.*V(2).*0.5./maxValue,...
%                     histHdl.BinEdges(i+1).*Vc(2)+Lp.*V(2).*1.05+histHdl.Values(i).*Lp.*V(2).*0.5./maxValue,...
%                     histHdl.BinEdges(i+1).*Vc(2)+Lp.*V(2).*1.05]./1.5,...
%                     [.8,0,0],'FaceAlpha',.5);
%             end

        else % less than 1500 per day will be excluded
            PCAslope9yrs(i)=nan;
        end
    end

figure;plot(1:3287,PCAslope9yrs);
% PCAslope=real(PCAslope);

disp('Calculated the PCA slope of each day of 9 years. Well Done!')

% save the file
save('D:\000learning\00research\00DataAndCode\Mengcheng_wkt_Nhindley\mchspPCA9yrs_2014_2022.mat', ...
    'PCAslope9yrs');
disp('Data 9yrs Saved! Please Check!')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the PCA slope of each day of ave year
% time use: about 5 mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0)
PCAslopeave=zeros(366,1);
    for i = 1:366
        id = find(doy >=i-1 & doy <i+1);
        if length(id)>1500 % more than 1500 per day will be included
            mch_id=mch(id);      speed_id=speed(id)./1000;
            % limit the data range and let
            % the length of X(:,1) and X(:,2) be the same
            fd= find(mch_id> 70 & mch_id < 110 ...
                & speed_id >8 & speed_id< 80);
            mch_id=mch_id(fd);
            speed_id=speed_id(fd);

            % X is the original data for PCA method
            X=zeros(200000,2);
            X(1:length(fd),1)=speed_id;
            X(1:length(fd),2)=mch_id;
%             n=1;
%             for sp=8:4:80
%                 m=1;
%                 for ht=70:2:110
%                     numdis_id(n,m)=length(find(mch_id>ht-1 & ...
%                         mch_id<=ht+1 & speed_id>sp-2 & speed_id<sp+2));
%                     m=m+1;
%                 end
%                 n=n+1;
%             end
%             
%             % numdis_id=numdis_id'
%             numdis_id=log(numdis_id/max(max(numdis_id)));
%             H1=scatter(X(:,1),X(:,2),'x','LineWidth',2,varargin{:},'visible','off');

% Calcu Main direction
            Xmean=mean(X);
            X0=X-Xmean;
            covMat=(X0.')*X0./size(X,1);
            [V,~]=eigs(covMat,1);
            V=V.*sign(V(1));
            Vc=[V(2);-V(1)];
            Lp=max(X0*V);
            Ln=min(X0*V);
            Ll=max(X0*Vc);
            Lr=min(X0*Vc);
% Draw Histogram
            H2=plot(Xmean(1)+[Lp,Ln].*V(1).*1.05,Xmean(2)+[Lp,Ln].*V(2).*1.05,'visible','off');

%% calcu the slope
            PCAslopeave(i)=(max(H2.YData)-min(H2.YData))/(max(H2.XData)-min(H2.XData));
%             % this is the plot of main direction
%             % hold on; plot(H2.XData, H2.YData)

%             % hold on;
%             % H3=fill(Xmean(1)+[Lp,Lp,Ln,Ln].*V(1).*1.05+[Ll,Lr,Lr,Ll].*Vc(1).*1.05,...
%             %      Xmean(2)+[Lp,Lp,Ln,Ln].*V(2).*1.05+[Ll,Lr,Lr,Ll].*Vc(2).*1.05,... 
%             %     H1.CData,'FaceAlpha',.15,'EdgeColor','none');
%             
%             histHdl=histogram(X0*Vc,'BinEdges',linspace(Lr,Ll,13),'Visible','off');
%             maxValue=max(histHdl.Values);
% 
%             for i=1:histHdl.NumBins
%                 H4(i)=fill(Xmean(1)+[histHdl.BinEdges(i).*Vc(1)+Lp.*V(1).*1.05,...
%                     histHdl.BinEdges(i).*Vc(1)+Lp.*V(1).*1.05+histHdl.Values(i).*Lp.*V(1).*0.5./maxValue,...
%                     histHdl.BinEdges(i+1).*Vc(1)+Lp.*V(1).*1.05+histHdl.Values(i).*Lp.*V(1).*0.5./maxValue,...
%                     histHdl.BinEdges(i+1).*Vc(1)+Lp.*V(1).*1.05]./1.5,...
%                     Xmean(2)+[histHdl.BinEdges(i).*Vc(2)+Lp.*V(2).*1.05,...
%                     histHdl.BinEdges(i).*Vc(2)+Lp.*V(2).*1.05+histHdl.Values(i).*Lp.*V(2).*0.5./maxValue,...
%                     histHdl.BinEdges(i+1).*Vc(2)+Lp.*V(2).*1.05+histHdl.Values(i).*Lp.*V(2).*0.5./maxValue,...
%                     histHdl.BinEdges(i+1).*Vc(2)+Lp.*V(2).*1.05]./1.5,...
%                     [.8,0,0],'FaceAlpha',.5);
%             end

        else % less than 1500 per day will be excluded
            PCAslopeave(i)=nan;
        end
    end

figure;plot(1:366,PCAslopeave);
% PCAslope=real(PCAslope);

disp('Calculated the PCA slope of each day of ave. Well Done!')

% save the file
save('D:\000learning\00research\00DataAndCode\Mengcheng_wkt_Nhindley\mchspPCAave_2014_2022.mat', ...
    'PCAslopeave');
disp('Data ave Saved! Please Check!')

end