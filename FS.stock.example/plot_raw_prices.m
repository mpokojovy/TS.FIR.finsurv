%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A Fast Initial Response Approach to Real-Time Financial Surveillance  %
%            (C) Michael Pokojovy and Andrews T. Anum (2022)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_raw_prices(dates, prices, ref_date, plot_title)
    figure;

    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 26; ySize = 12;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position', [0 0 xSize*50 ySize*50]);

    color_order = get(gca, 'colororder');
    %%
    subplot_tight(1, 1, 1, [0.08 0.05]);

    plot(dates, prices, '-', 'LineWidth', 2, 'Color', color_order(1, :));
    xline(ref_date);

    title(plot_title, 'interpreter', 'latex', 'FontSize', 18);

    xlabel('Date', 'interpreter', 'latex', 'FontSize', 18);
    ylabel('Closing price', 'interpreter', 'latex', 'FontSize', 18);
end