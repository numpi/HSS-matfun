function  H = full_to_hss2(A,options)
%Compute the standard hss decomposition of a full matrix A

if nargin<2
    options=[];
end

if isfield(options, 'tol')==0
    options.tol=1e-12;
end

if isfield (options, 'threshold')==0
    options.threshold=256;
end

if isfield (options, 'issymm')==0
    options.issymm=0;
end
if isfield (options, 'nrm')==0
    options.nrm=1;
end



m = options.threshold;

n=size(A,1);
if n < m
    H=A;
else
    H = hss2;
    H.n=n;
    H.leaf = 1;
    H.top=0;
    H.D = cell(1,ceil(n/m));
    H.U = cell(1,ceil(n/m));
    H.V = cell(1,ceil(n/m));
    for i = 1 : floor(n/m)
        H.D{i}=A((i-1)*m+1:i*m,(i-1)*m+1:i*m);
        [U,S]=svd(A((i-1)*m+1:i*m,[1:(i-1)*m,i*m+1:size(A,2)]));
        s=sum(diag(S)>options.tol*options.nrm*S(1,1));
        H.U{i}=U(:,1:s);
        if options.issymm
            H.D{i}=(H.D{i}+H.D{i}')/2;
            H.V{i}=H.U{i};
        else
        [U,S]=svd(A([1:(i-1)*m,i*m+1:size(A,2)],(i-1)*m+1:i*m)');
        s=sum(diag(S)>options.tol*options.nrm*S(1,1));
        H.V{i}=U(:,1:s);
        end
    end
    if floor(n/m)<ceil(n/m)
        i=ceil(n/m);
        H.D{i}=A((i-1)*m+1:end,(i-1)*m+1:end);
        [U,S]=svd(A((i-1)*m+1:end,1:(i-1)*m));
        s=sum(diag(S)>options.tol*options.nrm*S(1,1));
        H.U{i}=U(:,1:s);
        if options.issymm
            H.D{i}=(H.D{i}+H.D{i}')/2;
            H.V{i}=H.U{i};
        else
        [U,S]=svd(A(1:(i-1)*m,(i-1)*m+1:end)');
        s=sum(diag(S)>options.tol*options.nrm*S(1,1));
        H.V{i}=U(:,1:s);
        end
    end
    H = full_to_hss2_in(H,blkdiag(H.U{:})'*(A-blkdiag(H.D{:}))*blkdiag(H.V{:}),options);
end
end
%----------------------------------------

function  H = full_to_hss2_in(H,A,options)

l=length(H.U);

if l==2
    H.B=A;
    H.top=1;
else
    H.B = hss2;
    H.B.n=size(A,1);
    H.B.top=0;
    H.B.leaf = 0;
    H.B.D = cell(1,ceil(l/2));
    H.B.U = cell(1,ceil(l/2));
    H.B.V = cell(1,ceil(l/2));
    k2=0;
    for i=1:ceil(l/2)
        k1=k2+1;
        if i==ceil(l/2)
            k2=size(A,1);
        else
            k2=k1-1+size(H.U{1+2*(i-1)},2)+size(H.U{2+2*(i-1)},2);
        end
        H.B.D{i}=A(k1:k2,k1:k2);
        [U,S]=svd(A(k1:k2,[1:k1-1,k2+1:size(A,2)]));
        s=sum(diag(S)>options.tol*options.nrm*S(1,1));
        H.B.U{i}=U(:,1:s);
        if options.issymm
            H.B.D{i}=(H.B.D{i}+H.B.D{i}')/2;
            H.B.V{i}=H.B.U{i};
        else
        [U,S]=svd(A([1:k1-1,k2+1:size(A,1)],k1:k2)');
        s=sum(diag(S)>options.tol*options.nrm*S(1,1));
        H.B.V{i}=U(:,1:s);
        end
    end
    H.B = full_to_hss2_in(H.B,blkdiag(H.B.U{:})'*(A-blkdiag(H.B.D{:}))*blkdiag(H.B.V{:}),options);
end
end