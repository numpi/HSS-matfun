fprintf("\n size & timeCKR & timeCKM & timeDense & errCKR & errCKM \\\\ \n \\hline \n");
addpath('util')
rng(1);
m=50;
for n = [1024,2048,4096]
    [timeCKR,timeCKM,errCKR,errCKM,timeDense]=Comparison_inv_sqrt(n);
    fprintf('%d & %.2f & %.2f & %.2f & %1.2e& %1.2e \\\\ \n', n, timeCKR, timeCKM, timeDense, errCKR, errCKM);
end

for n = [8192,16384]
    [timeCKR,timeCKM]=Comparison_inv_sqrt(n,m);
    fprintf('%d & %.2f & %.2f & & &  \\\\ \n', n, timeCKR, timeCKM);
end

%fprintf("\n num. poles & timeCKR & timeCKM & errCKR & errCKM \\\\ \n \\hline \n");

%n=4096;
%for m = [10,20,40,80,160]
%    [timeCKR,timeCKM,errCKR,errCKM,timeDense]=Comparison_inv_sqrt(n,m);

%    fprintf('%d & %.2f & %.2f & %1.2e& %1.2e  \\\\ \n', m, timeCKR, timeCKM, errCKR, errCKM);
%end


