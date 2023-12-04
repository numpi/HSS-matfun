fprintf("\n size & blkrank & time {\\CKR} & time {\\CKM} & time {\\tt LM} & time ULV & time Dense & err \\CKR & err \\CKM & err {\\tt LM} & err ULV \\\\ \n \\hline \n");
addpath('util')
for n = [1024,2048,4096, 8192]
    [blkrank, timeCKR,timeCKM, timeML, timeInvHss,errCKR,errCKM, errML,errInvHss, timeDense]=Comparison_inv_banded(n);
    fprintf('%d & %d & %.2f & %.2f & %.2f & %.2f & %.2f & %1.2e & %1.2e & %1.2e &%1.2e\\\\ \n', n,blkrank, timeCKR, timeCKM, timeML, timeInvHss, timeDense, errCKR, errCKM, errML,errInvHss);
end

