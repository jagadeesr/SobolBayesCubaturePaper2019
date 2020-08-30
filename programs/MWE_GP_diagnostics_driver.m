gail.InitializeWorkspaceDisplay

%% Exponential Cosine example
fwh = 1;
dim = 3;
npts = 2^6;
nRep = 20;
nPlot = 2;
[~, rOptAll, thOptAll, fName] = ...
   MWE_gaussian_diagnostics_engine(fwh,dim,npts, ...
   [],[],nRep,nPlot);

%% Plot Exponential Cosine example
figH = figure;
plot(rOptAll, thOptAll, ...
   '.','MarkerSize',20,'color',MATLABBlue)
%axis([4 6 0.1 10])
%set(gca,'yscale','log')
title(['\(d = ' num2str(dim) ',\ n = ' num2str(npts) '\)'])
xlabel('Inferred \(r\)')
ylabel('Inferred \(\theta\)')
%print('-depsc',[fName '-rthInfer-n-' int2str(npts) '-d-' ...
%   int2str(dim)])
saveas(figH, [fName '-rthInfer-n-' int2str(npts) '-d-' ...
   int2str(dim) '.jpg'])

%% Tests with random function
rArray = [1.5 2 4];
nrArr = size(rArray,2);
fParArray = [0.5 1 2; 1 1 1; 1 1 1];
nfPArr = size(fParArray,2);
fwh = 3;
dim = 2;
npts = 2^6;
nRep = 20;
nPlot = 2;
thetaAll(nrArr,nfPArr) = 0;
rOptAll(nrArr,nfPArr,nRep) = 0;
thOptAll(nrArr,nfPArr,nRep) = 0;
for jjj = 1:nrArr
   for kkk = 1:nfPArr
      [thetaAll(jjj,kkk), rOptAll(jjj,kkk,:), thOptAll(jjj,kkk,:), fName] = ...
         MWE_gaussian_diagnostics_engine(fwh,dim,npts, ...
         rArray(jjj),fParArray(:,kkk),nRep,nPlot);
   end
end

%% Plot figures for random function
figH = figure
colorArray = {MATLABBlue,MATLABOrange,MATLABGreen,MATLABCyan,MATLABMaroon,MATLABPurple};
nColArray = length(colorArray);
for jjj = 1:nrArr
   for kkk = 1:nfPArr
      clrInd = mod(nfPArr*(jjj-1)+kkk-1,nColArray)+1;
      clr = colorArray{clrInd};
      plot(reshape(rOptAll(jjj,kkk,:),[nRep,1]), ...
         reshape(thOptAll(jjj,kkk,:),[nRep,1]), ...
         '.','MarkerSize',20,'color',clr)
      hold on
      scatter(rArray(jjj),thetaAll(jjj,kkk),200,clr,'s','filled')
   end
end
axis([1 6 0.01 100])
set(gca,'yscale','log')
title(['\(d = ' num2str(dim) ',\ n = ' num2str(npts) '\)'])
xlabel('Inferred \(r\)')
ylabel('Inferred \(\theta\)')
saveas(figH,[fName '-rthInfer-n-' int2str(npts) '-d-' ...
   int2str(dim) '.jpg'])
 
%% Keister example
fwh = 2;
dim = 3;
npts = 2^6;
nRep = 20;
nPlot = 2;
[~, rOptAll, thOptAll, fName] = ...
   MWE_gaussian_diagnostics_engine(fwh,dim,npts, ...
   [],[],nRep,nPlot);

%% Plot Keister example
figH = figure();
plot(rOptAll, thOptAll, ...
   '.','MarkerSize',20,'color',MATLABBlue)
%axis([4 6 0.5 1.5])
%set(gca,'yscale','log')
xlabel('Inferred \(r\)')
ylabel('Inferred \(\theta\)')
title(['\(d = ' num2str(dim) ',\ n = ' num2str(npts) '\)'])
saveas(figH,[fName '-rthInfer-n-' int2str(npts) '-d-' ...
   int2str(dim) '.jpg'])


%% Keister example
fwh = 2;
dim = 3;
npts = 2^10;
nRep = 20;
nPlot = 2;
[~, rOptAll, thOptAll, fName] = ...
   MWE_gaussian_diagnostics_engine(fwh,dim,npts, ...
   [],[],nRep,nPlot);

%% Plot Keister example
figH=figure();
plot(rOptAll, thOptAll, ...
   '.','MarkerSize',20,'color',MATLABBlue)
%axis([4 6 0.5 1.5])
%set(gca,'yscale','log')
xlabel('Inferred \(r\)')
ylabel('Inferred \(\theta\)')
title(['\(d = ' num2str(dim) ',\ n = ' num2str(npts) '\)'])
saveas(figH,[fName '-rthInfer-n-' int2str(npts) '-d-' ...
   int2str(dim) '.jpg'])




