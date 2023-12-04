function [H, time] = rational_evaluation(A,f,Z)
%Computes the matrix function r(A), where r is a rational approximation of
%f on Z.
[r, pol, res, zer, zj, fj, wj] = aaa(f, Z);
%err(end)
%length(pol)
%[res,pol] = weights(pol,f,Z);
n=size(A,1);
H=(f(1)-sum(res./(1-pol)))*hss('eye',n);
%H2=hss('zeros',n,n);

if nargout>1
    t=tic;
end

% for i=1:length(zj)
%     X=inv(A-zj(i)*hss('eye',n));
%     H=H+wj(i)*fj(i)*X;
%     H=compress(H);
%     H2=H2+wj(i)*X;
%     H2=compress(H);
% end 
%H=H/H2;

for i=1:length(pol)
    H=H+res(i)*inv(A-pol(i)*hss('eye',n));
    H=compress(H);
end 
if nargout>1
    time=toc(t);
end
end

%-------------------------------------

function [alpha,beta] = weights(beta,f,t)
%given a set of poles beta a function f and a set of points t in which we 
% want to approximate f by a rational function, the procedure computes the
% weights alpha, such that f(x) ~ sum_i alpha_i/(x-beta_i)

beta=beta(:);
t=t(:);
B=1./(t-beta.');
alpha=pinv(B)*f(t);
%alpha=B\f(t);
beta=beta(alpha~=0);
alpha=alpha(alpha~=0);
end

