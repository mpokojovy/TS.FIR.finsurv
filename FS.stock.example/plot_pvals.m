%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A Fast Initial Response Approach to Real-Time Financial Surveillance  %
%            (C) Michael Pokojovy and Andrews T. Anum (2022)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_pvals(dfs, pvals)
    figure;

    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 26; ySize = 12;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position', [0 0 xSize*50 ySize*50]);
    
    hold on;
    plot(dfs, pvals, 'LineWidth', 2);
    [~, ind] = max(pvals);
    plot(dfs(ind), pvals(ind), 'k*', 'LineWidth', 2, 'MarkerSize', 8); 
    hold off;
    
    title('One-Sample Kolmogorov-Smirnov Test for Student''s $t$ Distribution', 'interpreter', 'latex', 'FontSize', 18);

    xlabel('Degrees of freedom $\nu$', 'interpreter', 'latex', 'FontSize', 18);
    ylabel('$p$-value', 'interpreter', 'latex', 'FontSize', 18);
end