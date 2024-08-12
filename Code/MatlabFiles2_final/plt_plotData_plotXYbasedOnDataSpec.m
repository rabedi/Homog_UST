function plt_plotData_plotXYbasedOnDataSpec(inputPlotDataProp, x, y)

lw = inputPlotDataProp.val_lineWidth;
if (inputPlotDataProp.val_lineStyle == ':')
    lw = lw * 2;
end
if (strcmp(inputPlotDataProp.val_marker, 'none') == 1)
    plt1=plot(x, y, 'LineWidth', lw, 'Color', inputPlotDataProp.val_lineColor, 'LineStyle', inputPlotDataProp.val_lineStyle);
    plt1.Color(4)=inputPlotDataProp.val_transparency ;
else
    plot(x, y, 'LineWidth', lw, 'Color', inputPlotDataProp.val_lineColor, 'LineStyle', inputPlotDataProp.val_lineStyle, ...
        'Marker', inputPlotDataProp.val_marker, 'MarkerSize', inputPlotDataProp.val_markerSize, ...
        'MarkerFaceColor', inputPlotDataProp.val_markerFaceColor, 'MarkerEdgeColor', inputPlotDataProp.val_markerEdgeColor);
end
