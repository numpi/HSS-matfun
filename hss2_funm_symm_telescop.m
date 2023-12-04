function H = hss2_funm_symm_telescop(A,f,poles,options)
%H = hss2_funm_symm_telescop(A,f,poles)
%
% Computes a telescopic decomposition of an HSS approximation of f(A) where 
% A is in hss2 format
%
% Inputs    
% A                     is a matrix in hss2 format;
% f                     is handle functions
% poles                 are poles for the computations of f(A)
% options.threshold     is the threshold size under which f(A) is directly
%                       computed by diagonalization of A (default 256)
% options.isreal        if the matrix is real and the poles are symmetric 
%                       wrt the real axis, in such a case the algorithm 
%                       only use real arithmetic (default 'false')
% options.deflationTol  tolerance to deflate vectors in the computation of
%                       Krylov subspaces (default 1e-12)
% options.nofun         does not compute the D factors of the telescopic
%                       factorization (default 'false')

if nargin<4
    options=[];
end

if isfield (options, 'threshold')==0
    options.threshold=256;
end
if isfield (options, 'isreal')==0
    options.isreal=0;
end
if isfield (options, 'deflationTol')==0
    options.deflationTol =1e-12;
end

if isfield (options, 'nofun')==0
    options.nofun = false;
end

H=hss2;
H.leaf=A.leaf;
H.top=A.top;
H.n=A.n;

if mod(length(A.D),2)==1
    A.U{end+1}=[];
    A.V{end+1}=[];
    A.D{end+1}=[];
end

H.U=cell(1,length(A.U));
H.D=cell(1,length(A.D));
H.V=cell(1,length(A.U));
if~A.top
    A.B.n=0;
end
for i = 1:length(A.D)
    if size(A.D{i})==0
        H.D{i}=[];
        H.U{i}=[];
        H.V{i}=[];
    else
        if (length(poles))*size(A.U{i},2)>=size(A.D{i},1)
            H.U{i}=eye(size(A.D{i}));
            H.V{i}=eye(size(A.D{i}));
            H.D{i}=zeros(size(A.D{i}));
        else
            if abs(poles(1)) == inf
                [w,~] = qr(A.U{i},0);
                j=2;
                H.U{i} = w;
            elseif options.isreal && poles(2)==conj(poles(1))
                [w,~] = qr(A.U{i},0);
                w = (A.D{i} - poles(1)*eye(size(A.D{i}))) \  w;
                w1=real(w);
                w=imag(w);
                j=3;
                [H.U{i},~] = qr([w1, w],0);
                w=H.U{i}(:,end-size(A.U{i},2)+1:end);
            else
                [w,~] = qr(A.U{i},0);
                w = (A.D{i} - poles(1)*eye(size(A.D{i}))) \  w;
                j=2;
                [w,~] = qr(w,0);
                H.U{i} = w;
            end

            while j <= length(poles)
                if abs(poles(j)) == inf
                    w=A.D{i} * w;
                    j=j+1;
                elseif options.isreal && j+1<=length(poles) && poles(j+1)==conj(poles(j))
                    w = (A.D{i} - poles(j)*eye(size(A.D{i}))) \  w;
                    w1=real(w);
                    w=imag(w);
                    w1 = w1 - H.U{i} * H.U{i}'*w1;
                    [w1,~] = qr(w1,0);
                    [H.U{i},~] = qr([H.U{i}, w1],0);
                    j=j+2;
                else
                    w = (A.D{i} - poles(j)*eye(size(A.D{i}))) \  w;
                    j=j+1;
                end
                [H.U{i},S]=svd([H.U{i}, w],0);
                s=sum(diag(S)>options.deflationTol *S(1,1)); 
                H.U{i}=H.U{i}(:,1:s);
                w=H.U{i}(:,end-size(A.U{i},2)+1:end);
            end
            H.V{i}=H.U{i};

            D=A.D{i};
            
            if options.nofun
                H.D{i}=A.D{i};
            else
            H.D{i}=funm(D,f)-H.U{i}*(funm(H.U{i}'*D*H.V{i},f))*H.V{i}';
            H.D{i}=(H.D{i}+H.D{i}')/2;
            end
        end
    end
    if mod(i,2)==0
        if A.top
            A.B=blkdiag(H.U{:})'*(blkdiag(A.U{:})*A.B*blkdiag(A.V{:})'+blkdiag(A.D{:}))*blkdiag(H.V{:});
            A.B=(A.B+A.B')/2;
        else
            A.B.D{i/2}=blkdiag(H.U{i-1:i})'*(blkdiag(A.U{i-1:i})*A.B.D{i/2}*blkdiag(A.V{i-1:i})'+blkdiag(A.D{i-1:i}))*blkdiag(H.V{i-1:i});
            A.B.D{i/2}=(A.B.D{i/2}+A.B.D{i/2}')/2;
            A.B.n=A.B.n+size(A.B.D{i/2},1);
            A.B.U{i/2}=blkdiag(H.U{i-1:i})'*blkdiag(A.V{i-1:i})*A.B.V{i/2};
            A.B.V{i/2} = A.B.U{i/2};
        end
    end
end


H.V=H.U;
if A.top
    if options.nofun
     H.B=A.B;
            else
    H.B=funm(A.B,f);
    H.B=(H.B+H.B')/2;
    end
elseif A.B.n<=options.threshold
    if options.nofun
                H.B=A.B;
            else
    H.B=funm(full(A.B),f);
    H.B=(H.B+H.B')/2;
    H = full_to_hss2_in(H,H.B,options);
    end
else
    A.B.V=A.B.U;
    [H.B] = hss2_funm_symm_telescop(A.B,f,poles,options);
end
end



%--------------------------------
function F=funm(A,f)
[V, D] = eig(A);
D = diag(f(diag(D)));
F = V * D / V;
end

function A=stab(A,U)
V=null(U');
A11=U'*A*U;
A22=V'*A*V;
A12=U'*A*V;
A21=V'*A*U;
S=A11-A12*(A22\A21);
m=min(eig(S));
A=A-(m-1)*(U*U');
end

function  H = full_to_hss2_in(H,A,options)


if isfield(options, 'deflationTol')==0
    options.deflationTol =1e-12;
end

if isfield (options, 'threshold')==0
    options.threshold=256;
end

options.issymm=1;

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
        s=sum(diag(S)>options.deflationTol *S(1,1));
        H.B.U{i}=U(:,1:s);
        if options.issymm
            H.B.V{i}=H.B.U{i};
        else
            [U,S]=svd(A([1:k1-1,k2+1:size(A,2)],k1:k2)');
            s=sum(diag(S)>options.deflationTol *S(1,1));
            H.B.V{i}=U(:,1:s);
        end
    end
    H.B = full_to_hss2_in(H.B,blkdiag(H.B.U{:})'*(A-blkdiag(H.B.D{:}))*blkdiag(H.B.V{:}),options);
end
end
