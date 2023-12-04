fprintf("n & err \\CKR & err {\\tt expm} &time \\CKR &time {\\tt expm} \\\\ \n \\hline \n");
addpath('util')
rng(1);
for n=2.^[10:13]

h=1/(n-1);
A=-spdiags(ones(n,1) * [-1 2 -1], -1:1, n, n);
A=(1/h^2)*A;


opts.r=2;
opts.m=4;
opts.treshold=256;
H = full_to_hss2(full(A),opts);
eps=1e-8;
hssoption('threshold',eps);
options.deflationTol =1e-8;


   
FA=expm(full(A));
nrm=norm(FA,"fro");

    poles=poles_exp2;
    t1=tic;
    FH1 = hss2_funm_symm_telescop(H,@(x) exp(x),[inf, poles],opts);
    timeCKR=toc(t1);

    errCKR=norm(full(FH1)-FA,"fro")/nrm;

    hssA=hss(A);
    t3=tic;
    FH3 = expm(hssA);
    timeexpm=toc(t3);
    errexpm=norm(full(FH3)-FA,"fro")/nrm;


    %fprintf('$%d$ & %1.2e & %1.2e & %.2f & %.2f\\\\ \n', n,errCKR,errCKM, timeCKR , timeCKM);
    fprintf('$%d$ & %1.2e & %1.2e & %.2f & %.2f\\\\ \n', n,errCKR,errexpm, timeCKR , timeexpm);
end

%-------------------------------

function poles = poles_exp2
v=[1,7.2270585462136953087e-1,2.5829131995630046565e-1,6.0809507390072969545e-2,...
    1.0598023489381217368e-2,1.4566258090943959608e-3,1.6421910387421514909e-4,...
    1.5592669068294270482e-5, 1.2701896118955897646e-6, 8.9952033533504590178e-8...
    5.5866285320057224635e-9, 3.0690410653109898959e-10, 1.4763544085912781412e-11,...
    6.5692146982166977577e-13, 2.1929711332561360156e-14,9.9803140529347369232e-16...
    1.0076443559489818370e-17, 9.9533336444995010637e-19 ];
m=length(v);

% no balancing
%C=diag(ones(1,m-2),-1)+diag(ones(1,m-2),1);
%C(m-1,m-2)=2;
%C(1,:)=-v(end-1:-1:1)./v(end);

%balancing
 C=8*diag(ones(1,m-2),-1)+diag(ones(1,m-2),1)/8;
 C(m-1,m-2)=2*8;
 C(1,:)=-v(end-1:-1:1)./v(end); 
 for i=1:size(C,1)
      C(1,i)=8^(-i+1)*C(1,i);
 end

C=C./2;
poles=-eig(C).';
end
