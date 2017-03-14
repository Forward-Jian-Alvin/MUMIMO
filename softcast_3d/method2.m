function [ chunkNum ] = method2( Lamb1,Lamb2,Lamb3,r)
%METHOD2 Summary of this function goes here
%   Detailed explanation goes here
    SLa1 = sum(sum(sum(Lamb1)));
    SLa2 = sum(sum(sum(Lamb2)));
    SLa3 = sum(sum(sum(Lamb3)));
    Lamb1 = Lamb1./SLa1;
    Lamb2 = Lamb2./SLa2;
    Lamb3 = Lamb3./SLa3;
    
    temp1 = reshape(Lamb1,1,[]);
    temp1 = sort(temp1,'descend');
    temp2 = reshape(Lamb2,1,[]);
    temp2 = sort(temp2,'descend');
    temp3 = reshape(Lamb3,1,[]);
    temp3 = sort(temp3,'descend');
 
    temp = [temp1 temp2 temp3];
    
    [tempp Index1] = sort(temp,'descend'); 
    Nn = numel(temp)*r;
    Index = Index1(1:Nn);
    
    ind1 = find(Index<=64*8);
    ind2 = find(Index>64*8 & Index<=64*8*2);
    ind3 = find(Index>64*8*2);
    
    chunkNum = [numel(ind1) numel(ind2) numel(ind3)];

end

