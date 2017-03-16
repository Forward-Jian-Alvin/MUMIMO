function [ a ] = IDCT3( A )
%IDCT3 Summary of this function goes here
%   Detailed explanation goes here

[H W Z] = size(A);
a = zeros(size(A));

for i = 1 : Z
    a(:,:,i) = idct(idct(A(:,:,i)')');
end

if Z ~= 1
a = reshape(a, H*W, Z);
a = idct(a');
a = a';
a = reshape(a, H, W, Z);
end
% for i = 1 : W
%     for j = 1 : H
%         t = idct(a(i,j,:));
%         a(i,j,:) = t;
%     end
% end
end

