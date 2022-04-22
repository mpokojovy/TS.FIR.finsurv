%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A Fast Initial Response Approach to Real-Time Financial Surveillance  %
%            (C) Michael Pokojovy and Andrews T. Anum (2022)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pvals, df_best] = student_KS_test(x, dfs)
    n = length(x);
    pvals = zeros(size(dfs));
    
    for ind = 1:length(dfs)
        df = dfs(ind);

        H0_cdf = makedist('tlocationscale', 'mu', median(x), ... 
                          'sigma', (quantile(x, 0.75) - quantile(x, 0.25))/(tinv(0.75, df) - tinv(0.25, df)), 'nu', df);
        [~, pval] = kstest(sort(x), 'CDF', H0_cdf);
        pvals(ind) = pval;
    end
    
    [~, i_best] = max(pvals);
    df_best = dfs(i_best);
end