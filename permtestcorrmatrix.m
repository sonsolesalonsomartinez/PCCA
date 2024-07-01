function pval = permtestcorrmatrix(Cin,Nperm)
% tests if the diagonal of the correlation matrix is higher than the off-diagonal
% Diego Vidaurre (2021)

    N = length(Cin);
    if (nargin<2) || isempty(Nperm)
        Nperm = 10000;
    end

    grotperms = zeros(Nperm,1);
    Cin0 = Cin;
    for perm=1:Nperm
        if (perm>1)
            Cin = Cin0(randperm(N),:);
        end
        grotperms(perm) = sum(diag(Cin));
    end

    pval = sum(grotperms>=grotperms(1)) / (Nperm+1);

end
