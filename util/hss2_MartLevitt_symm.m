function H = hss2_MartLevitt_symm(A,options)
%H = hss2_MartLevitt_symm(A,r,m,options)
%
% Computes a telescopic decomposition of a Hermitian HSS matrix A as 
% described in [1]. The algorithm is based only on matrix-vector products 
% v-->A*v where v is a random (block) vector. 
%
% Inputs
% A                 can be either a matrix or a handle function that performs
%                   the matrix-vector product v-->A*v, in the latter case 
%                   options.sizeA is needed;
%
%options.tol        required tolerance (default 1e-12)
%options.m          is the starting threshold size of A. In needed, this 
%                   value is adaptively increased by the algorithm (default 15)
%
%options.r          is the starting block rank size of A. In needed, this 
%                   value is adaptively increased by the algorithm
%                   (default floor(options.m/2))
%options.oversampl  is the oversampling parameter (defaulf 5)
%
%options.sizeA      number of columns of A. Needed only if A is a handle
%                   function
%
% Output
% H                 hss2 representation of A

n = size(A,1);

if nargin<2
    options=[];
end

if isfield(options, 'tol')==0
    options.tol=1e-12;
end

if isfield (options, 'oversampl')==0
    options.oversampl=5;
end

if isfield (options, 'm')==0
    options.m=15;
end

if isfield (options, 'r')==0
    options.r=floor(options.m/2);
end

if ~isa(A,'function_handle')
    n=size(A);
    A=@(v) A*v;
else 
    if isfield (options, 'sizeA')==0
        error('The number of columns of A has to be provided in options.sizeA')
    else
        n=options.sizeA;
    end
end

%-----------------------------------------

check=false;


r=options.r;
m=options.m;
Omega = randn(n,m+r+options.oversampl);
Y=A(Omega);

while ~check
    T=cell(1,ceil(n/m));
    D=cell(1,ceil(n/m));

    for i = 1:floor(n/m)
        P=null(Omega((i-1)*m+1:i*m,:));
        P=P(:,1:r+options.oversampl);
        [T{i},S] = svd(Y((i-1)*m+1:i*m,:)*P);
        T{i}=T{i}(:,1:sum(abs(diag(S))>options.tol*S(1,1)));
        if size(T{i},2)>r
            check=false;
            r=size(T{i},2);
        else
            check= true;
        end
        if ~check
            break
        else
            r=max(r,size(T{i},2));
        end

        D{i}=(eye(size(T{i},1))-T{i}*T{i}')*(Y((i-1)*m+1:i*m,:)/Omega((i-1)*m+1:i*m,:));
        D{i}=D{i}+T{i}*T{i}'*D{i}';
    end
    if check
        i=length(T);
        P=null(Omega((i-1)*m+1:end,:));
        P=P(:,1:min(r+options.oversampl,size(P,2)));
         [T{i},S] = svd(Y((i-1)*m+1:end,:)*P);
        T{i}=T{i}(:,1:sum(abs(diag(S))>options.tol*S(1,1)));

        if size(T{i},2)>r
            check=false;
            r=size(T{i},2);
        else
            check= true;
        end

        if check

        D{i}=(eye(size(T{i},1))-T{i}*T{i}')*Y((i-1)*m+1:end,:)/Omega((i-1)*m+1:end,:);
        D{i}=D{i}+T{i}*T{i}'*D{i}';
        end
    end
    if ~check
        m=2*r;
        k=3*r;
        Omega2=randn(n,k-size(Omega,2)+options.oversampl);
        Omega=[Omega,Omega2];
        Y=[Y,A(Omega2)];
    else
        if m<2*r
            m=2*r;
            k=3*r;
            Omega2=randn(n,k-size(Omega,2)+options.oversampl);
            Omega=[Omega,Omega2];
            Y=[Y,A (Omega2)];
            check=false;
        else
            H=hss2;
            H.U=T;
            H.V=T;
            H.D=D;
            H.leaf=1;
            H.top=0;
            [H,check,r] = hss2_small_MartLevitt_symm(H,r,Omega,Y,options);
            if ~check
                m=2*r;
                k=3*r;
                Omega2=randn(n,k-size(Omega,2)+options.oversampl);
                Omega=[Omega,Omega2];
                Y=[Y,A(Omega2)];
            end
        end
    end
end
%disp(['required matvec = ', num2str(size(Omega,2))])
end

%-----------------------------------------------

function [H,check,r] = hss2_small_MartLevitt_symm(H,r,Omega,Y,options)

if isfield(options, 'tol')==0
    options.tol=1e-12;
end

Y=blkdiag(H.U{:})'*(Y-blkdiag(H.D{:})*Omega);
Omega=blkdiag(H.V{:})'*Omega;
l=length(H.U);
if l==2
    H.B=Y/Omega;
    H.B=(H.B+H.B')/2;
    H.top=1;
    check=true;
else
    if mod(l,2)==1
        H.U{l+1}=[];
        H.V{l+1}=[];
        H.D{l+1}=[];
    end
    H.B=hss2;
    H.B.n=size(Omega,1);

    T=cell(1,ceil(l/2));
    D=cell(1,ceil(l/2));
    k2=0;
    for i=1:ceil(l/2)
        k1=k2+1;
        if i==ceil(l/2)
            k2=size(Omega,1);
        else
            k2=k1-1+size(H.U{1+2*(i-1)},2)+size(H.U{2+2*(i-1)},2);
        end

        P=null(Omega(k1:k2,:));
        P=P(:,1:min(r+options.oversampl,size(P,2)));
        [T{i},S] = svd(Y(k1:k2,:)*P);
        T{i}=T{i}(:,1:sum(abs(diag(S))>options.tol*S(1,1)));
        
            if size(T{i},2)>r
                check=false;
                r=size(T{i},2);
            else
                check= true;
            end
            if ~check
                H=0;
                r=max(size(T{i},2),r+1);
                return
            end

        D{i}=(eye(size(T{i},1))-T{i}*T{i}')*(Y(k1:k2,:)/Omega(k1:k2,:));
        D{i}=D{i}+T{i}*T{i}'*D{i}';
    end

    H.B.U = T;
    H.B.V = T;
    H.B.D = D;
    H.B.leaf = 0;
    H.B.top = 0;

    [H.B,check,r] = hss2_small_MartLevitt_symm(H.B,r,Omega,Y,options);

end
end






