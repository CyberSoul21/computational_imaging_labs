%load HDR
% hdr = hdrread('../Part1/Results/hdr_image.hdr');
hdr = hdrread('../Part1/Results/ours_hdr_image.hdr');

fignum = 6;

for dR = [2, 4, 6]

    % imwrite( ...
    %     durand_tonemapping(hdr, dR, fignum), ...
    %     "Results/tonemapped_durand_dR_" + dR + ".png" ...
    %     )
    imwrite( ...
        durand_tonemapping(hdr, dR, fignum), ...
        "Results/ours_tonemapped_durand_dR_" + dR + ".png" ...
        )
    
    fignum = fignum + 1;

end

clear fignum dR;
