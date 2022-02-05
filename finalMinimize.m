AbsorptionData = xlsread('/Users/liguanda/Desktop/MasterReasearch/LGD_Journal_Paper/NEW/AllSpectrum.xlsx');
AbsorptionData = AbsorptionData(104:end,1:end);
freq = AbsorptionData(:,1);
AbsorptionData = AbsorptionData(:,2:end);
[totalInformation,txtInfo,~] = xlsread('/Users/liguanda/Desktop/MasterReasearch/LGD_Journal_Paper/NEW/AllInformation.xlsx');

Onethird = [125 160 200,250,315,400,500,630,800,1000,1250,1600];
Indx = zeros(1,length(Onethird));
for n=1:1:length(Onethird)
    Indx(1,n) = find (freq == Onethird(1,n));
end

x0 = [1.5,5e3,0.5];
lb = [1,1e3,0];
ub = [3,8e4,5];

TestSample = length(AbsorptionData(1,:));
DataPoolPade = cell(7,TestSample);

for num = 1:TestSample
    DataPoolPade{1,num} = txtInfo{1,num};
    DataPoolPade{2,num} = AbsorptionData(:,num);
    
    rng default
    % options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
    % without Pore
    % fun = @(x)JCAmodel(x(1),x(2),Information(9,num),x(3),x(4),Information(7,num),freq)-TotalAbsorptionData(:,(num+1));
    
    % with Pore
    fun = @(x)H&S_model(x(1),x(2),totalInformation(9,num),x(3),totalInformation(7,num),freq)-AbsorptionData(:,(num));
    
    problem = createOptimProblem('lsqnonlin','x0',x0,'objective',fun,'lb',lb,'ub',ub);
    ms = MultiStart('PlotFcns',@gsplotbestf);
    [x,errormultinonlin] = run(ms,problem,100);
    % 
    % [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub);
    
    
    DataPoolPade{3,num} = H&S_model(x(1),x(2),totalInformation(9,num),x(3),totalInformation(7,num),freq);
    
    DataPoolPade{4,num} = x;
    DataPoolPade{5,num} = x(1,1);
    DataPoolPade{6,num} = x(1,2);
    DataPoolPade{7,num} = DataPoolPade{3,num}-DataPoolPade{2,num};
    DataPoolPade{8,num} = 1-sqrt(sum(DataPoolPade{7,num}.^2)/sum(DataPoolPade{2,num}.^2));    
    
    
    figure(num);
    plot(freq,DataPoolPade{2,num});
    hold on;
    plot(freq,DataPoolPade{3,num});
    hold off;
    legend('test','predict');
end

%% plot tort to stress
figure(101)
tortS = zeros(8,5);
tortS_x = [0 0.1 0.2 0.3 0.4];
for j = 1:1:5
    tortS(1,j) = DataPoolPade{5,4*j-3};
    tortS(2,j) = DataPoolPade{5,4*j-2};
    tortS(3,j) = DataPoolPade{5,4*j-1};
    tortS(4,j) = DataPoolPade{5,4*j};
    tortS(5,j) = sum(tortS(1:4,j))/4;
    tortS(6,j) = max(tortS(1:4,j))-tortS(5,j);
    tortS(7,j) = tortS(5,j)-min(tortS(1:4,j));
end 
plot(tortS_x,tortS(5,:));
errorbar(tortS_x,tortS(5,:),tortS(7,:),tortS(6,:),'-o','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
ylim([1 3]);
xlim([0 0.4]);
ylabel('曲折度因子(-)','fontsize',12);
xlabel('成型压强(MPa)','fontsize',12);
%set(gcf,'position',[200,200,480,360]);
grid on
%% plot tort to agg0.2
figure(102)
tortA2 = zeros(8,5);
tortA2_x = [0.7 0.8 0.9 1.0 1.1];
for k = 1:1:5
    tortA2(1,k) = DataPoolPade{5,4*k+17};
    tortA2(2,k) = DataPoolPade{5,4*k+18};
    tortA2(3,k) = DataPoolPade{5,4*k+19};
    tortA2(4,k) = DataPoolPade{5,4*k+20};
    tortA2(5,k) = sum(tortA2(1:4,k))/4;
    tortA2(6,k) = max(tortA2(1:4,k))-tortA2(5,k);
    tortA2(7,k) = tortA2(5,k)-min(tortA2(1:4,k));
end 
plot(tortA2_x,tortA2(5,:));
errorbar(tortA2_x,tortA2(5,:),tortA2(7,:),tortA2(6,:),'-o','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
ylim([1 3]);
xlim([0.7 1.1]);
ylabel('曲折度因子(-)','fontsize',12);
xlabel('骨胶比','fontsize',12);
%set(gcf,'position',[200,200,480,360]);
grid on

%% plot tort to agg0.4
figure(103)
tortA4 = zeros(8,5);
tortA4_x = [0.7 0.8 0.9 1.0 1.1];
for m = 1:1:5
    tortA4(1,m) = DataPoolPade{5,4*m+37};
    tortA4(2,m) = DataPoolPade{5,4*m+38};
    tortA4(3,m) = DataPoolPade{5,4*m+39};
    tortA4(4,m) = DataPoolPade{5,4*m+40};
    tortA4(5,m) = sum(tortA4(1:4,m))/4;
    tortA4(6,m) = max(tortA4(1:4,m))-tortA4(5,m);
    tortA4(7,m) = tortA4(5,m)-min(tortA4(1:4,m));
end 
plot(tortA4_x,tortA4(5,:));
errorbar(tortA4_x,tortA4(5,:),tortA4(7,:),tortA4(6,:),'-o','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
ylim([1 3]);
xlim([0.7 1.1]);
ylabel('曲折度因子(-)','fontsize',12);
xlabel('骨胶比','fontsize',12);
%set(gcf,'position',[200,200,480,360]);
grid on
%% plot sigma to stress 
figure(201)
sigmaS = zeros(8,5);
sigmaS_x = [0 0.1 0.2 0.3 0.4];
for j = 1:1:5
    sigmaS(1,j) = DataPoolPade{6,4*j-3};
    sigmaS(2,j) = DataPoolPade{6,4*j-2};
    sigmaS(3,j) = DataPoolPade{6,4*j-1};
    sigmaS(4,j) = DataPoolPade{6,4*j};
    sigmaS(5,j) = sum(sigmaS(1:4,j))/4;
    sigmaS(6,j) = max(sigmaS(1:4,j))-sigmaS(5,j);
    sigmaS(7,j) = sigmaS(5,j)-min(sigmaS(1:4,j));
end 
plot(sigmaS_x,sigmaS(5,:));
errorbar(sigmaS_x,sigmaS(5,:),sigmaS(7,:),sigmaS(6,:),'-o','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
%ylim([1 3]);
xlim([0 0.4]);
ylabel('流阻(N?s/m^2)','fontsize',12);
xlabel('成型压强(MPa)','fontsize',12);
%set(gcf,'position',[200,200,480,360]);
grid on
%% plot sigma to agg0.2
figure(202)
sigmaA2 = zeros(8,5);
sigmaA2_x = [0.7 0.8 0.9 1.0 1.1];
for j = 1:1:5
    sigmaA2(1,j) = DataPoolPade{6,4*j+17};
    sigmaA2(2,j) = DataPoolPade{6,4*j+18};
    sigmaA2(3,j) = DataPoolPade{6,4*j+19};
    sigmaA2(4,j) = DataPoolPade{6,4*j+20};
    sigmaA2(5,j) = sum(sigmaA2(1:4,j))/4;
    sigmaA2(6,j) = max(sigmaA2(1:4,j))-sigmaA2(5,j);
    sigmaA2(7,j) = sigmaA2(5,j)-min(sigmaA2(1:4,j));
end 
plot(sigmaA2_x,sigmaA2(5,:));
errorbar(sigmaA2_x,sigmaA2(5,:),sigmaA2(7,:),sigmaA2(6,:),'-o','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
%ylim([1 3]);
xlim([0.7 1.1]);
ylabel('流阻(N?s/m^2)','fontsize',12);
xlabel('骨胶比(-)','fontsize',12);
%set(gcf,'position',[200,200,480,360]);
grid on
%% plot sigma to agg0.4
figure(203)
sigmaA4 = zeros(8,5);
sigmaA4_x = [0.7 0.8 0.9 1.0 1.1];
for j = 1:1:5
    sigmaA4(1,j) = DataPoolPade{6,4*j+37};
    sigmaA4(2,j) = DataPoolPade{6,4*j+38};
    sigmaA4(3,j) = DataPoolPade{6,4*j+39};
    sigmaA4(4,j) = DataPoolPade{6,4*j+40};
    sigmaA4(5,j) = sum(sigmaA4(1:4,j))/4;
    sigmaA4(6,j) = max(sigmaA4(1:4,j))-sigmaA4(5,j);
    sigmaA4(7,j) = sigmaA4(5,j)-min(sigmaA4(1:4,j));
end 
plot(sigmaA4_x,sigmaA4(5,:));
errorbar(sigmaA4_x,sigmaA4(5,:),sigmaA4(7,:),sigmaA4(6,:),'-o','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
%ylim([1 3]);
xlim([0.7 1.1]);
ylabel('流阻(N?s/m^2)','fontsize',12);
xlabel('骨胶比(-)','fontsize',12);
%set(gcf,'position',[200,200,480,360]);
grid on
%%
figure(301)
a21 = plot(freq,DataPoolPade{2,33},'--o','MarkerIndices',Indx,'MarkerSize',3,'Color',[0.3 0.3 0.3]);
set(a21,'MarkerFaceColor',get(a21,'color'));
hold on;
a22 = plot(freq,DataPoolPade{3,33},'-o','MarkerIndices',Indx,'MarkerSize',3,'Color',[0 0 0]);
set(a22,'MarkerFaceColor',get(a22,'color'));
hold on;
a23 = plot(freq,DataPoolPade{2,35},'--o','MarkerIndices',Indx,'MarkerSize',3,'Color',[1 0 0]);
set(a23,'MarkerFaceColor',get(a23,'color'));
hold on;
a24 = plot(freq,DataPoolPade{3,35},'-o','MarkerIndices',Indx,'MarkerSize',3,'Color',[0.6350 0.0780 0.1840]);
set(a24,'MarkerFaceColor',get(a24,'color'));
hold off;
LocalLegend = legend('试件1实测','试件1拟合','试件2实测','试件2拟合','Location','Northwest');
xlim([100 1600]);
ylim([0 1]);
ylabel('吸声系数 \alpha','fontsize',12);
xlabel('频率(Hz)','fontsize',12);
set(gcf,'position',[200,200,480,360]);
set(gca,'XScale','log');
set(LocalLegend,'box','off');
grid on;
set(gca,'xtick',Onethird);
