function [L, DroppedN] = fx_DropLeastE64Coefficient(Cdct, ToDropN)

[H, W, GOP] = size(Cdct);

blkH = H / 8;
blkW = W / 8;
chunksize = blkH * blkW;

Chunk = zeros(64, chunksize, GOP);
lambda = zeros(GOP, 64);
for indFrm = 1:GOP
    Chunk(:,:,indFrm) = fx_Reshape_Img2Block(Cdct(:,:,indFrm), blkH, blkW);
    lambda(indFrm, :) = mean(Chunk(:,:,indFrm).^2, 2);
end

DroppedN = 0;
while DroppedN < ToDropN
    [~, indDropChunk] = min(lambda(:));
    [indFrm, indChunk] = ind2sub(size(lambda), indDropChunk);
    
    if ToDropN - DroppedN > chunksize
        DropN = chunksize;
    else
        DropN = ToDropN - DroppedN;
    end
    Chunk(indChunk, end-DropN+1:end, indFrm) = -inf;
    DroppedN = DroppedN + DropN;
    lambda(indDropChunk) = inf;
end

LinChunk = zeros(size(Chunk));
for indFrm = 1:GOP
    for indChunk = 1:64
        Coeff = Chunk(indChunk, :, indFrm);
        mask = find(Coeff~=-inf);
        if isempty(mask)
            LinChunk(indChunk, :, indFrm) = 0;
        else
            LinChunk(indChunk, mask, indFrm) = mean(Coeff(mask).^2);
        end
    end
end

L = zeros(size(Cdct));
for indFrm = 1:GOP
    L(:,:,indFrm) = fx_Reshape_Block2Img(LinChunk(:,:,indFrm), blkH, blkW, H, W);
end

end