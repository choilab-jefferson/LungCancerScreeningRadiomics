function [ax, cor, sag] = extract_three_slices(I, slices)
    ax = I(:,:,slices(1));
    cor = squeeze(I(:,slices(2),:));
    sag = squeeze(I(slices(3),:,:));
end