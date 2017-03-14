function [ VarVid_Y ] = VarEqualChunk( Y, BlockNum )
    fun4 = @(x)  sum((x(:)-0).^2)/length(x(:)) * ones(size(x));
    VarVid_Y = zeros(size(Y));
    [W, H, Z] = size(Y);
    BandBlkListInfo = getMultiLevelDWTInfoEqualSizeBlk_pro(W, H, BlockNum);
    for j = 1:Z
        VarVid_Y(:,:,j)  = BandInplaceProc(Y(:,:,j), BandBlkListInfo, fun4);
    end
end



