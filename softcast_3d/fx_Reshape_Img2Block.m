function block = fx_Reshape_Img2Block(img, block_height, block_width)

    [height, width, GOP] = size(img);

    MapImg2Blk = fx_CreateImgBlockMap(height, width, block_height, block_width);
    [Nb, bsz] = size(MapImg2Blk);
    block = zeros(Nb * GOP, bsz);
    for ii = 1:GOP
        tmp = img(:,:,ii);
        block((ii-1)*Nb+(1:Nb), :) = tmp(MapImg2Blk);
    end
end