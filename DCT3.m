function [ A ] = DCT3( a )
%DCT3 Summary of this function goes here
%   Detailed explanation goes here
[H, W, Z] = size(a);
A = zeros(size(a));

for i = 1 : Z
    A(:,:,i) = dct(dct(a(:,:,i))')';
end

if Z ~= 1
    A = reshape(A, H*W, Z);
    A = dct(A');
    A = A';
    A = reshape(A, H, W, Z);
end

% for i = 1 : W
%     for j = 1 : H
%         A(i,j,:) = dct(A(i,j,:));
%     end
% end

end