function color_array = ColorRGB_ForPlot(color_indx)

% This is a small function which saves 6 different colour schemes used for
% the R2* paper (and other papers). The first three are for in silico data
% at variable g-ratio, while the next three are for ex vivo data at
% variable fibre dispersion. The last three are only to the colour bars in
% Figure 8 (\epsilon_m). The colour schemes were obtained from colorbrew.

    switch color_indx
        case 1 % In silico g-ratio 0.66
            color_array = [158 202 225]/256;
        case 2 % In silico g-ratio 0.73
            color_array = [49 130 189]/256;
        case 3 % In silico g-ratio 0.8
            color_array = [8 81 156]/256;
        case 4 % Ex vivo highly dispersed
            color_array = [161 217 155]/256;
        case 5 % Ex vivo mildly dispersed
            color_array = [49 163 84]/256;
        case 6 % Ex vivo negligible dispersed
            color_array = [0 109 44]/256;
        case 7
            color_array = [217 217 217]/256;
        case 8
            color_array = [150 150 150]/256;
        case 9
            color_array = [99 99 99]/256;
    end
end