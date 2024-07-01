function [pval,grotperms,r,U,V,stats,A,B] = permtestcca0(Xin,Yin,NPCA,Nperm,confounds,cs,Permutations,ps)
% tests if the canonical correlation between Xin and Yin is different from zero
% Diego Vidaurre, University of Oxford (2021)

N = length(Yin);
if nargin<3, NPCA = []; end
if (nargin>4) && ~isempty(confounds)
    confounds = confounds - repmat(mean(confounds),N,1);
    Xin = Xin - confounds * pinv(confounds) * Xin;
    Yin = Yin - confounds * pinv(confounds) * Yin;
end
if (nargin<8)
    ps = []; %permutation structure
end
if (nargin<6)
    cs = []; %family structure
elseif ~isempty(cs)
    [grotMZi(:,2),grotMZi(:,1)]=ind2sub([length(cs) length(cs)],find(tril(cs,1)==1));
    [grotDZi(:,2),grotDZi(:,1)]=ind2sub([length(cs) length(cs)],find(tril(cs,1)==2));
end
if (nargin<4) || isempty(Nperm)
    Nperm = 10000;
end
if (nargin<7) || isempty(Permutations)
    Permutations = [];
end
PrePerms=0;
if ~isempty(Permutations)
    PrePerms=1;
    Nperm=size(Permutations,2);
end

Xin = Xin - mean(Xin); Xin = Xin ./ std(Xin);
Yin = Yin - mean(Yin); Yin = Yin ./ std(Yin);
if ~isempty(NPCA)
    if length(NPCA) == 1, NPCA = [NPCA NPCA]; end
    if NPCA(1) > 0
        if any(isnan(Xin(:)))
            [~,Xin] = pca(Xin, 'Algorithm','als');
        else
            [~,Xin] = pca(Xin);
        end
        Xin = Xin(:,1:NPCA(1));
    end
    if NPCA(2) > 0
        if any(isnan(Yin(:)))
            [~,Yin] = pca(Yin, 'Algorithm','als');
        else
            [~,Yin] = pca(Yin);
        end
        Yin = Yin(:,1:NPCA(2));
    end
end

grotperms = zeros(Nperm,1); Yin0 = Yin;
for perm=1:Nperm
    if (perm>1)
        if PrePerms==1 % pre-supplied permutation
            Yin=Yin0(Permutations(:,perm),:);  % or maybe it should be the other way round.....?
        elseif isempty(cs) && isempty(ps)      % simple full permutation with no correlation structure
            rperm = randperm(N);
            Yin=Yin0(rperm,:);
        elseif isempty(cs) && ~isempty(ps)     % within-subject permutation with no correlation structure
            PERM=zeros(1,N);
            subID = unique(ps);
            for j = 1:length(subID)
                ind = find(ps == subID(j));
                rp = randperm(length(ind));
                PERM(ind) = ind(rp);
            end
            Yin=Yin0(PERM,:);
        elseif ~isempty(cs)                     % 3.between-subject permutation with correlation structure
            PERM=zeros(1,length(cs));
            perm1=randperm(size(grotMZi,1));
            for ipe=1:length(perm1) % permute within MonoZig families
                if rand<0.5, wt=[1 2]; else wt=[2 1]; end
                PERM(grotMZi(ipe,1))=grotMZi(perm1(ipe),wt(1));
                PERM(grotMZi(ipe,2))=grotMZi(perm1(ipe),wt(2));
            end
            perm1=randperm(size(grotDZi,1)); % permute within DiZig families
            for ipe=1:length(perm1)
                if rand<0.5, wt=[1 2]; else wt=[2 1]; end
                PERM(grotDZi(ipe,1))=grotDZi(perm1(ipe),wt(1));
                PERM(grotDZi(ipe,2))=grotDZi(perm1(ipe),wt(2));
            end
            from=find(PERM==0);  pto=randperm(length(from));  to=from(pto);  PERM(from)=to;
            if isempty(ps)                      % 3.a. simple
                Yin=Yin0(PERM,:);
            else                                % 3.b. for all trials
                perm1=PERM;PERM=zeros(1,N);
                for j = 1:length(perm1)
                    ind = ps == perm1(j);
                    PERM(ind) = find(ps == j);
                end
                Yin=Yin0(PERM,:);
            end
        end
    end

    [A,B,r,U,V,stats] = canoncorr(Xin,Yin);
    grotperms(perm) = sum(abs(r));
end

if any(isnan(grotperms(:))), error('NaN appeared..'); end
pval = sum(grotperms>=grotperms(1)) / (Nperm+1);

end