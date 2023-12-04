fprintf("\n $n$ & time \\CKR & time \\CKM & time {\\tt LM} & time InvHss & time Dense & err \\CKR & err \\CKM & err {\\tt LM} & err {\\tt inv}\\\\ \n \\hline \n");
addpath('util')
for n = [1024,2048,4096, 8192]
    [timeCKR,timeCKM, timeML, timeInvHss ,errCKR,errCKM, errML,errInvHss, timeDense]=Comparison_inv(n);
    fprintf('%d & %.2f & %.2f & %.2f & %.2f & %.2f & %1.2e & %1.2e & %1.2e & %1.2e \\\\ \n', n, timeCKR, timeCKM, timeML, timeInvHss, timeDense, errCKR, errCKM, errML, errInvHss);
end

for n = [16384, 32768]
    [timeCKR,timeCKM,timeML, timeInvHss]=Comparison_inv(n);
    fprintf('%d & %.2f & %.2f & %.2f & %.2f & & & & &  \\\\ \n', n, timeCKR, timeCKM,timeML,timeInvHss);
end
