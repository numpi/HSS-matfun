function [blkrank,timeCKR,timeCKM, timeML, timeInvHss ,errCKR,errCKM, errML,errInvHss, timeDense]=Comparison_inv_banded(n)
% Computation of sqrtm(-A)^(-1)  where A is a HSS matrix

rng(1);
A = GL_matrix(n, 1.5);
blkrank= hssrank(A);

opts.treshold = 256;
opts.tol = 1e-14;
opts.issymm = true;

H = full_to_hss2(full(A),opts);
%H = hss2_to_standard(H);
f=@(x) 1./x;
t1=tic;
FH=hss2_funm_symm_telescop(H,f,0);
timeCKR=toc(t1);

t2=tic;
  FH2 =  CKM2(A, @(x) inv(x), 0);
  timeCKM=toc(t2);

opts.sizeA=n;
opts.m=4 * hssrank(A);
opts.r=2 * hssrank(A);
t3=tic;
F=ulv(A);
AA=@(v) ulv_solve(F,v);
FH3 = hss2_MartLevitt_symm(AA,opts);
timeML=toc(t3);

t4=tic;
FH4 = inv(A);
timeInvHss=toc(t4);

if nargout > 5
t4=tic;
FA=inv(full(A));
timeDense=toc(t4);
nrm=norm(full(FA),"fro");
errCKR=norm(full(FA)-full(FH),"fro")/nrm;
errCKM=norm(full(FA)-full(FH2),"fro")/nrm;
errML=norm(full(FA)-full(FH3),"fro")/nrm;
errInvHss=norm(full(FA)-full(FH4),"fro")/nrm;
end
end

%--------------------------------------------

function A = GL_matrix(n, alpha)
%GL_matrix n is the size, 1 <= alpha <= 2 the fractional power.

[am, ap] = my_fractional_symbol(alpha, n);
ap(n) = 0;
A = hss('toeplitz', am, ap);
A = A + A';
A = symmetrize(A, 'up');

end

function [am, ap] = my_fractional_symbol(alpha, n)
%FRACTIONAL_SYMBOL Construct the symbol of the Grunwald-Letkinov derivative
%
% [AM, AP] = FRACTIONAL_SYMBOL(ALPHA, N) construct the negative and
%     positive parts of the symbol of the Toeplitz matrix discretizing the
%     fractional derivative by means of the Grunwald-Letkinov shifted
%     formulas. 
%

v = zeros(n+2, 1);
v(1) = 1;

v = -cumprod([1,1-((alpha + 1)./(1:n))]);
am = v(2:end);
ap = [v(2), v(1)];

end
