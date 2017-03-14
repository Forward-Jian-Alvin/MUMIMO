function [Info]=getMultiLevelDWTInfoEqualSizeBlk_pro(H, W, BlockNum)

n = sqrt(BlockNum);
assert(mod(H, n)  == 0);
assert(mod(W, n)  == 0);

H_blk = H / n;
W_blk = W / n;

fields = {'BlockIdx', 'Top', 'Left', 'Height', 'Width'};
c = cell(BlockNum, length(fields));

BlockIdx = 0;
for i=1:n
    for j=1:n
        c{BlockIdx+1,1} = BlockIdx+1;
        c{BlockIdx+1,2} = (i-1) * H_blk + 1;
        c{BlockIdx+1,3} = (j-1) * W_blk + 1;
        c{BlockIdx+1,4} = H_blk;
        c{BlockIdx+1,5} = W_blk;
        BlockIdx = BlockIdx +1;
    end
end
Info = cell2struct(c, fields, 2);