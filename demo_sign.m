n=4096;
addpath('util')
rng(1);
fprintf("\n a & time \\CKR & time \\CKM & time Rational &  err \\CKR & err \\CKM & errRational \\\\ \n \\hline \n");
for a=-1:-2:-9
b=0;
eps=1e-12;
l= logspace(a,b,n/2);
l=[-l(end:-1:1),l];
X=randn(n);
[X,~]=qr(X);
A=X*diag(l)*X';
[P,A] = hess(A);
A=tril(A)+tril(A,-1)';

FA = P'*X*diag(sign(l))*X'*P;

nrm=norm(FA,"fro");

opts.issymm=1;
opts.isreal=1;
opts.threshold=256;
H = full_to_hss2(A);

f=@(x) sign(real(x));
poles = poles_Zolot_sign(10^a,10^b,eps);
t=tic;
FH1 = hss2_funm_symm_telescop(H,f,[inf,poles],opts);
timeCKR=toc(t);

errCKR=norm(full(FH1)-FA,"fro")/nrm;

HssA=hss(A);
t=tic;
FH6 = CKM2(hss(A),@(x)signm(x),poles);
timeCKM=toc(t);

errCKM=norm(full(FH6)-FA,"fro")/nrm;


f=@(x) sign(real(x));
[FH3,timeRat] = rational_evaluation(HssA,f,l);
errRational=norm(full(FH3)-FA,"fro")/nrm;


 fprintf('%d & %.2f & %.2f & %.2f & %1.2e& %1.2e & %1.2e  \\\\ \n', a, timeCKR, timeCKM,timeRat, errCKR, errCKM, errRational);
end

%-----------------------------------

function poles = poles_Zolot_sign(a,b,tol)
%Computes m optimal zolotarev poles for the sign function on [-b,-a]U[a,b]
mu=(1-sqrt(a/b))/(1+sqrt(a/b));
m=ceil(log(2/tol + 1) / (pi*ellipke(sqrt(1-mu^2))/(4*ellipke(mu))));
sq=sqrt(1-(a/b)^2);
K = ellipke(sq);
s=ceil(m/2);
poles=zeros(1,2*s);
for i=1:s
    sn=ellipj((2*i-1)*K/(2*s),sq);
    c=(a)*sqrt(-sn^2/(1-sn^2));
    poles(2*i-1)=c;
    poles(2*i)=-c;
end
end

function F = signm (A)
[V, D] = eig(A);
D = diag(sign(real((diag(D)))));
F = V * D / V;
end

