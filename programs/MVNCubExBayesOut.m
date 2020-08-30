function MVNCubExBayesOut(dataFileName)

% MVNCubExBayesOut: Prints LaTeX table of outputs from
% MVNCubatureExampleBayes

gail.InitializeDisplay

dirpath = '.';
load(dataFileName);
%% Output
fid = fopen([dirpath, filesep, 'MVNBayesOut.txt'], 'wt');

fprintf(fid,'\\[ \n \\begin{array}{rccccc} \n \\toprule \n');
istart = 2; 
for i = istart:length(avgAbsErrMC)
   if i == 1
      fprintf(fid,'\\multicolumn{6}{c}{\\varepsilon = \\num{%2.0e}, \\ d = %1.0f,\\ \\mSigma = \\mI, \\ \\vb=-\\va=(3.5,\\dots,3.5) } \\\\ \\midrule \n', abstol, d);
   else
	   fprintf(fid,'\\multicolumn{6}{c}{\\varepsilon = \\num{%2.0e}, \\ d = %1.0f,\\ \\mSigma = 0.4 \\mI, + 0.6\\vone\\vone^T , \\ \\va=(-\\infty,\\dots,-\\infty), \\ \\vb\\sim \\sqrt{d}\\cu[0,1]^d } \\\\ \\midrule \n',abstol,d);
   end
  fprintf(fid,' \\text{Method} & \\text{MC} & \\text{Lattice} & \\text{Sobol''} & \\text{BayesLat} & \\text{BayesSobol}  \\\\\n');
  fprintf(fid,' \\text{Absolute Error} & \\num{%3.2e} & \\num{%3.2e} & \\num{%3.2e}  & \\num{%3.2e}  & \\num{%3.2e}  \\\\\n', ...
     round([avgAbsErrMC(i); avgAbsErrLat(i); avgAbsErrSob(i); avgAbsErrBayLat(i); avgAbsErrBaySob(i)], 2, 'significant'));
  fprintf(fid,' \\text{Tolerance Met} & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%%  \\\\\n', ...
     100*round([succMC(i); succLat(i); succSob(i); succBayLat(i); succBaySob(i)], 2, 'significant'));
  fprintf(fid,' n & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} \\\\\n', ...
     round([nSampleMC(i); nSampleLat(i); nSampleSob(i); nSampleBayLat(i); nSampleBaySob(i)], 2, 'significant'));
  fprintf(fid,' \\text{Time (seconds)} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f}  \n', ...
     round([timeMC(i); timeLat(i); timeSob(i); timeBayLat(i); timeBaySob(i)], 2,'significant'));
  if i == 1, fprintf(fid,' \\\\ \\\\ \n'); end
end
fprintf(fid,'\\end{array} \n \\]');

fclose(fid);

end
