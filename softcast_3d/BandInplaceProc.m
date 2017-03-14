function [Result]=BandInplaceProc(X, BandListInfo, func)
Result = zeros(size(X));
nBands = length(BandListInfo);
for i=1:nBands
    RowList = BandListInfo(i).Top  : BandListInfo(i).Top  + BandListInfo(i).Height - 1;
    ColList = BandListInfo(i).Left : BandListInfo(i).Left + BandListInfo(i).Width  - 1;
    Result(RowList, ColList)=func(X(RowList, ColList));
end