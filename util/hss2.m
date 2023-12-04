classdef hss2
    %hss decomposition described in [1]
    %
    %[1] Levitt, J., & Martinsson, P. G. (2022). Linear-complexity black-box
    %randomized compression of hierarchically block separable matrices.
    %arXiv preprint arXiv:2205.02990.

    properties
        % The matrix is decomposed as
        %      blkdiag(U{:})*B*blkdiag(V{:})'+blkdiag(D{:})
        %
        %
        U
        V
        D
        B

        n %size of the matrix

        top
        leaf
    end
    methods
        function A=full(obj)
            if isequal(class(obj),'hss2')
                if obj.top
                    A=blkdiag(obj.U{:})*obj.B*blkdiag(obj.V{:})'+blkdiag(obj.D{:});
                else
                    A=blkdiag(obj.U{:})*full(obj.B)*blkdiag(obj.V{:})'+blkdiag(obj.D{:});
                end
            else
                A=obj;
            end
        end
    end
end