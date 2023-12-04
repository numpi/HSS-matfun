function [timeCKR,timeCKM, timeML, timeInvHss,errCKR,errCKM, errML,errInvHss, timeDense]=Comparison_inv(n)
% Computation of A^(-1) where A the discretization of the
%Laplacian of size n.

rng(1);

h=1/(n-1);
A=spdiags(ones(n,1) * [-1 2 -1], -1:1, n, n);
A=(1/h^2)*A;

opts.treshold=256;
opts.tol=1e-14;

H = full_to_hss2(full(A),opts);
H = hss2_to_standard(H);
f=@(x) 1./x;
t1=tic;
FH=hss2_funm_symm_telescop(H,f,0);
timeCKR=toc(t1);

HssA=hss(A);
t2=tic;
FH2 = CKM2(HssA,@(x) inv(x),0);
timeCKM=toc(t2);

opts.sizeA=n;
opts.m=4;
opts.r=2;
t3=tic;
F=ulv(HssA);
AA=@(v) ulv_solve(F,v);
FH3 = hss2_MartLevitt_symm(AA,opts);
timeML=toc(t3);

t4=tic;
FH4 = inv(HssA);
timeInvHss=toc(t4);

if nargout > 4
t4=tic;
FA=inv(A);
timeDense=toc(t4);
nrm=norm(full(FA));
errCKR=norm(FA-full(FH))/nrm;
errCKM=norm(FA-full(FH2))/nrm;
errML=norm(FA-full(FH3))/nrm;
errInvHss=norm(FA-full(FH4))/nrm;
end
end
