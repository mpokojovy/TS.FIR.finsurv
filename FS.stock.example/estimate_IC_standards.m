%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A Fast Initial Response Approach to Real-Time Financial Surveillance  %
%            (C) Michael Pokojovy and Andrews T. Anum (2022)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu0, sigma0] = estimate_IC_standards(x, est_type, par)
    est_type = lower(est_type);

    if (strcmp(est_type, 'mdpd'))
        alpha_mdpd = par;
        [~, mu0, sigma0] = ddiv_estimator(mean(x), std(x), x, alpha_mdpd, false);
    elseif (strcmp(est_type, 'mcd'))
        bdp = par;
        mcd = mcd1D(x, par);
        mu0    = mcd.loc;
        sigma0 = mcd.cov;
    else
        mu0    = mean(x);
        sigma0 = std(x);
    end
end