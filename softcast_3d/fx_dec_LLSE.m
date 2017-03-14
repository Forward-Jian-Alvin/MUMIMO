function C = fx_dec_LLSE(Y, G, L, Hc, sigma_n, hadamardsize)

if isempty(Hc)
    rho = G.*L;
    C = rho./(G.*rho + sigma_n^2) .* Y;
else
    GOP   = numel(Y);
    C     = cell(GOP, 1);
    blksz = size(Hc, 2);
    
    if hadamardsize == 64
        Ha = hadamard(64) / 8;
    else
        Ha = 1;
    end
    
    indDiag = find(eye(size(G{1}, 1)));
    
    for ii = 1:GOP
        h      = Hc(:,:,ii);
        y      = Y{ii};
        g      = G{ii};
        lambda = L{ii};
        Block = [];
        for jj = 1:blksz
            hc = diag(h(:,jj));
            hr = [real(hc) -imag(hc); imag(hc) real(hc)];
            MatA  = hr * Ha * diag(g(:,jj));
            La = diag(lambda(:,jj));

            yr = [real(y(:,jj)); imag(y(:,jj))];
            
            yr(yr == 1e100) = 0;

            %blk = La * A' * pinv( A * La * A' + N0 * eye(size(A,1))) * yr;
            MatB = La * MatA';
            MatC = MatA * MatB;
            MatC(indDiag) = MatC(indDiag) + sigma_n^2;
            blk = MatB / MatC * yr;

            Block = [Block blk];
        end
        C{ii} = Block;
    end
end