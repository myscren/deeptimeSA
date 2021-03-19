function varargout = LSA(varargin)
% This program is to implment the 1D local singularity analysis (LSA) using
% power-law scaling transform. The input dataset is a two dimensional array,
% data(:,1) represents the time/locations while data(:,2) for the value. 
% Begin initialization code
% clc;clear all;
%% input dataset
[filename1 filepath1] = uigetfile({'*.xlsx';'*.xls'},'Select grd file');
input_path = strcat(filepath1,filename1);
data = xlsread(input_path);
data2 = data(:,2);
% -min(data(:,2))+5; % plus one in order to avoid log(0) in Line 22

%% multiscale analysis using moving avarage
winds=[3:11]; % the range of the window size
for scale=1:length(winds)
    MAdata(scale,:)=smooth(data2,winds(scale));  
end

%% estimate local singularity index using scaling transform
for i=1:length(data2)
    for scale=1:length(winds)
        log_value(1,scale)=log(MAdata(scale,i)*winds(scale)); 
        log_rscale(1,scale)=log(winds(scale));
    end
    [a1,a2]=m2c(log_rscale,log_value); % least square fitting
    delta_alpha(i,1)=1-a1;
    fra_dens(i,1)=exp(a2);
   % calculate the fitting goodness
    f_fit=log_rscale.*a1+a2;
    cor=corrcoef(f_fit,log_value);
    R2(i,1)=cor(1,2); % the coefficient of determination
end

%% plot
figure(1)
set(gcf,'color','white');grid on
subplot(2,1,1)
plot(data(:,1),data(:,2));hold on;
plot(data(:,1),fra_dens,'r'); axis tight;
xlabel('Mys','FontSize',10);
ylabel('Numbers','FontSize',10);
legend('Original data','Fractal density');
subplot(2,1,2)
plot(data(:,1),delta_alpha,'r'); axis tight;grid on
xlabel('Mys','FontSize',110);
ylabel('Singularity','FontSize',10);

%% save result
result=[data(:,1),delta_alpha];
xlswrite('singularity_500_sed_area',result); % save singularity value
xlswrite('R2_500_sed_area',R2); % save the coefficient of determination


function [a1,a2]=m2c(x,y)
% m2c(x,y) performs the least square fitting for y=a1*x+a2 
% [x,y] is the input data, a1 and a2 are the estimates of slope and intercept, respectively.  
Sumx=0.0; 
Sumy=0.0;
Sumxy=0.0;
Sumx2=0.0;
a1=0.0;
a2=0.0;
m=length(x);    
for i=1:m 
	Sumx=Sumx+x(i)*1.0;
	Sumy=Sumy+y(i)*1.0;
	Sumxy=Sumxy+(x(i)*y(i))*1.0;
	Sumx2=Sumx2+(x(i)*x(i))*1.0;
end    
a1=(m*Sumxy-Sumx*Sumy)/(m*Sumx2-Sumx*Sumx);%%% slope
a2=(Sumx2*Sumy-Sumx*Sumxy)/(m*Sumx2-Sumx*Sumx);%%% intercept

%% references
% [1]Cheng Q. Singularity analysis of global zircon U-Pb age series and
% implication of continental crust evolution. Gondwana Research, 2017,
% 51:51-63. 
% [2]Cheng Q. Mapping singularities with stream sediment geochemical
% data for prediction of undiscovered mineral deposits in Gejiu, Yunnan
% Province, China. Ore Geology Reviews, 2007, 32(1–2):314-324.
% [3]Chen, G., and Q. Cheng, Fractal density modeling of crustal
% heterogeneity from the KTB deep hole, Journal of Geophysical Research
% Solid Earth, 2017, 122.
