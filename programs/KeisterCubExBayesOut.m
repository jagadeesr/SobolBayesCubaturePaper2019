% KeisterCubExBayesOut: Prints LaTeX table of outputs from
% KeisterCubatureExampleBayes
function KeisterCubExBayesOut(dataFileName)

gail.InitializeDisplay

dirpath = '.';
load(dataFileName);
istart =  1;
ifinish = 1;
%% Output
fid = fopen([dirpath, filesep, 'KeisterBayesOut.txt'], 'wt');
fprintf(fid,'\\[ \n \\begin{array}{rccccc} \n \\toprule \n');
for i = 1:ifinish
  fprintf(fid,'\\multicolumn{6}{c}{\\varepsilon = \\num{%2.0e}, \\ d = %1.0f} \\\\ \n', abstol(i), dvec(i));
  fprintf(fid,' \\hline \n');
  fprintf(fid,' \\text{Method} & \\text{MC} & \\text{Lattice} & \\text{Sobol} & \\text{BayesLat} & \\text{BayesSobol}  \\\\ \n');
  fprintf(fid,' \\text{Absolute Error} & \\num{%1.6e} & \\num{%1.2e} & \\num{%1.2e}  & \\num{%1.2e}  & \\num{%1.2e}  \\\\ \n', ...
     round([avgAbsErrMC(i); avgAbsErrLat(i); avgAbsErrSob(i); avgAbsErrBayLat(i); avgAbsErrBaySob(i)], 2, 'significant'));
  fprintf(fid,' \\text{Tolerance Met} & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%%  \\\\\n', ...
     100*round([succMC(i); succLat(i); succSob(i); succBayLat(i); succBaySob(i)], 2, 'significant'));
  fprintf(fid,' n & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} \\\\\n', ...
     round([nSampleMC(i); nSampleLat(i); nSampleSob(i); nSampleBayLat(i); nSampleBaySob(i)], 2, 'significant'));
  fprintf(fid,' \\text{Time (seconds)} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} \\\\ \n', ...
     round([timeMC(i); timeLat(i); timeSob(i); timeBayLat(i); timeBaySob(i)], 2,'significant'));
  fprintf(fid,' \\\\ \n');
end
fprintf(fid,'\\end{array} \n \\]');

fclose(fid);


