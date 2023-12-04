function [timeCKR,timeCKM,errCKR,errCKM,timeDense,timeML,errML,errBound]=Comparison_inv_sqrt(n,m)
% Computation of sqrtm(-A)^(-1) performing m steps of rational
%Krylov with quasi optimal poles where A the discretization of the
%Laplacian of size n.

rng(1);

h=1/(n-1);
A=spdiags(ones(n,1) * [-1 2 -1], -1:1, n, n);
A=(1/h^2)*A;
opts.treshold=256;

e=eig(A);
a=min(e); b=max(e);
opts.issymm=1;
opts.nrm=b;
H = full_to_hss2(full(A),opts);

if nargin<2
    eps=1e-12;
    Phi=@(x)(((2.*x-(b+a))/(b-a))-sqrt(((2.*x-(b+a))/(b-a)).^2-1));
    m=ceil(log(8*log2(n/opts.treshold)/(sqrt(a)*eps*abs(Phi(0))))*log(16*b/a)/pi^2);
end

f=@(x) 1./sqrt(x);
poles = poles_Markov_functions(a,b,-Inf,0,m);
opts.deflationTol=hssoption('threshold');
t1=tic;
FH=hss2_funm_symm_telescop(H,f,[inf,poles],opts);
timeCKR=toc(t1);
HssA=hss(A);

t2=tic;
FH2 = CKM2(HssA,@(x) inv(sqrtm(x)),poles);
timeCKM=toc(t2);

if nargout > 2
    t4=tic;
    FA=inv(sqrtm(full(H)));
    timeDense=toc(t4);
    nrm=norm(full(FA));
    errCKR=norm(FA-full(FH))/nrm;
    errCKM=norm(FA-full(FH2))/nrm;
end
if nargout>5
    opts.sizeA=n;
    AA=@(v) krylov_fAb(H2,v,f,poles);
    t3=tic;
    FH3 = hss2_MartLevitt_symm(AA,opts);
    timeML=toc(t3);
    errML=norm(FA-full(FH3))/nrm;
end
if nargout > 7
    errBound=16/sqrt(a)*expm(-m*pi^2/(log(16*b/a)));
end
end

%-----------------------------------------

function v = krylov_fAb(A,z,f,poles)
v=zeros(size(z));
for i=1:5:size(z,2)
    if size(z,2)<i+5
        b=z(:,i:end);
    else
        b=z(:,i:i+5);
    end
    [V,H,K]= rat_krylov(A, b, [poles,inf]);
    k=size(b,2);
    if cond(K(1:end-k, :)) < 1e10 % was 15
        Ap = H(1:end-k, :) / K(1:end-k, :);
    else
        warning('MartLevitt:: ill-conditioned matrix K, projection computed with additional matvecs')
        Ap = V(:, 1:end-k)' * A * V(:, 1:end-k);
    end
    Ap=(Ap+Ap')/2;
    [U,D] = eig(full(Ap));
    D = diag(f(diag(D)));
    F = U * D / U;
    if size(z,2)<i+5
        v(:,i:end) = V(:,1:end-k)*F*(V(:,1:end-k)'*b);
    else
        v(:,i:i+5) = V(:,1:end-k)*F*(V(:,1:end-k)'*b);
    end

end
end
%-----------------------------------------------
function poles=poles_Markov_functions(a,b,alpha,beta,m)
%eigenvalues in [a,b] integration domain (alpha,beta)
phi=@(x)(((2.*x-(b+a))/(b-a))-sqrt(((2.*x-(b+a))/(b-a)).^2-1));
if abs(alpha)==Inf
    k=-1/phi(beta);
else
    k=(phi(alpha)-phi(beta))/(1-phi(alpha)*phi(beta));
end
k=((1-sqrt(1-k^2))/k)^2;
%T1=@(x)(1-phi(beta).*x)./(x-phi(beta));
T1=@(x)(1+phi(beta)*x)/(x+phi(beta));
T2=@(x) (sqrt(k).*x-1)./(x-sqrt(k));
T=@(x) T1(T2(x));
K = ellipke(k);
poles=zeros(1,m);
for i=1:m
    poles(i)=sqrt(k)*ellipj(K*(m+1-2*i)/m,k);
    poles(i)=T(poles(i));
end
end

