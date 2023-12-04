function H= hss2_to_standard(HH)
%function H= hss2_to_standard(HH)
%
% converts any telescopic representation into the standard one
%
%Input HH  hss2
%Output H  hss2

H=HH;
if mod(length(H.D),2)==1
    H.D{end+1}=[];
    H.U{end+1}=[];
    H.V{end+1}=[];
end
if H.top==1
    for i=1:length(H.D)
        H.D{i}=find_diag_block(H,i);
    end
    i=2;
    H.B=HH.B-blkdiag(H.U{i-1:i})'*(-blkdiag(HH.D{i-1:i})+blkdiag(H.D{i-1:i}))*blkdiag(H.U{i-1:i});
else
    for i=1:length(H.D)
        H.D{i}=find_diag_block(HH,i);
        if mod(i,2)==0
            HH.B.D{i/2}=HH.B.D{i/2}-blkdiag(H.U{i-1:i})'*(-blkdiag(HH.D{i-1:i})+blkdiag(H.D{i-1:i}))*blkdiag(H.U{i-1:i});
        end
    end
    H.B=hss2_to_standard(HH.B);
end
end

%-------------------------------------------------------

function D = find_diag_block(A,i)
s=size(A.U{i},2);
if A.top
    if i==2        
        D=A.D{i}+A.U{i}*A.B(end-s+1:end,end-s+1:end)*A.U{i}';
    else
        D=A.D{i}+A.U{i}*A.B(1:s,1:s)*A.U{i}';
    end
else
D=A.D{i}+A.U{i}*(find_diag_block_in(A.B,i,s))*A.U{i}';
end
end

%---------------------------------------------------------------------

function D = find_diag_block_in(A,i,s)
b=mod(i,2);
i=ceil(i/2);
if A.top
    if b==0
        D=A.D{i}(end-s+1:end,end-s+1:end);
        if i==1
           D=D+A.U{i}(end-s+1:end,:)*A.B(1:size(A.U{i},2),1:size(A.U{i},2))*A.U{i}(end-s+1:end,:)';
        else
           D=D+A.U{i}(end-s+1:end,:)*A.B(end-size(A.U{i},2)+1:end,end-size(A.U{i},2)+1:end)*A.U{i}(end-s+1:end,:)';
        end
    else
        D=A.D{i}(1:s,1:s);
        if i==1
           D=D+A.U{i}(1:s,:)*A.B(1:size(A.U{i},2),1:size(A.U{i},2))*A.U{i}(1:s,:)';
        else
           D=D+A.U{i}(1:s,:)*A.B(end-size(A.U{i},2)+1:end,end-size(A.U{i},2)+1:end)*A.U{i}(1:s,:)';
        end
    end
else
    if b==0
        D=A.D{i}(end-s+1:end,end-s+1:end);
        D=D+A.U{i}(end-s+1:end,:)*(find_diag_block_in(A.B,i,size(A.U{i},2)))*A.U{i}(end-s+1:end,:)';
    else
        D=A.D{i}(1:s,1:s);
        D=D+A.U{i}(1:s,:)*(find_diag_block_in(A.B,i,size(A.U{i},2)))*A.U{i}(1:s,:)';
    end
end
end
