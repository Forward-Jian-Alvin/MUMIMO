function img = fx_Reshape_Block2Img(block, block_height, block_width, height, width, GOP)

img = zeros(height, width, GOP);

MapImg2Blk = fx_CreateImgBlockMap(height, width, block_height, block_width);

[Nb, bsz] = size(MapImg2Blk);

for ii = 1:GOP
    tmp = zeros(height, width);
    tmp(MapImg2Blk) = block((ii-1)*Nb+(1:Nb), :);
    img(:,:,ii) = tmp;
end
