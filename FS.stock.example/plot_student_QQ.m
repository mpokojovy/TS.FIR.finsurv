%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A Fast Initial Response Approach to Real-Time Financial Surveillance  %
%            (C) Michael Pokojovy and Andrews T. Anum (2022)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_student_QQ(x, dfs, nrow, ncol)
    if (length(dfs) ~= nrow*ncol)
        error("length(dfs) must be equal nrow*ncol");
    end

    figure;
    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 28; ySize = 18;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position', [0 0 xSize*50 ySize*50]);

    n = length(x);
    pvals = zeros(size(dfs));

    for ind = 1:length(dfs)
        df = dfs(ind);

        H0_cdf = makedist('tlocationscale', 'mu', median(x), ... 
                          'sigma', (quantile(x, 0.75) - quantile(x, 0.25))/(tinv(0.75, df) - tinv(0.25, df)), 'nu', df);
        [~, pval] = kstest(sort(x), 'CDF', H0_cdf);
        pvals(ind) = pval;

        subplot_tight(nrow, ncol, ind, [0.1 0.05]);

        hold on;
        plot(tinv(((1:n) - 0.5)/n, df), sort(x), 'b+');
        xlim([min(tinv(((1:n) - 0.5)/n, df)), max(tinv(((1:n) - 0.5)/n, df))]);
        ylim(0.5*(min(x) + max(x)) + 0.75*(max(x) - min(x))*[-1 1]);

        xgrid = linspace(tinv(1E-4, df), tinv(1 - 1E-4, df), 2);
        ygrid = median(x) + (quantile(x, 0.75) - quantile(x, 0.25))/(tinv(0.75, df) - tinv(0.25, df))*xgrid;

        plot(xgrid, ygrid, 'r--');

        title(['Robust Student''s $t_{', num2str(df), '}$ Q-Q-plot'], 'interpreter', 'latex', 'FontSize', 18);

        xlabel('Theoretical quantiles', 'interpreter', 'latex', 'FontSize', 18);
        ylabel('Empirical quantiles', 'interpreter', 'latex', 'FontSize', 18);
    end
end