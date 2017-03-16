function fx_MIMOchannelInit(Ngop, Nc)

    global mimoChan;

    fname = sprintf('meta_MIMOfading_%d_%d.mat', Ngop, Nc);
    if exist(fname, 'file') == 2
        load(fname, 'H');
    else
        H = 1/sqrt(4) * (randn (4, 4, Ngop, Nc) + 1i*randn (4, 4, Ngop, Nc));
%        H = cell(Ngop, Nc);
%        for jj = 1:Ngop
%            for ii = 1:Nc
%                H{jj,ii} = 1/sqrt(2) * (randn (2, 2) + 1i*randn (2, 2));
%            end
%        end
        save(fname, 'H');
    end

    if size(H, 3) < Ngop || size(H, 4) < Nc
        error('fx_MIMOchannelInit: less channel fading gains for simulation.');
    end

    mimoChan.Hrayleigh = H;
    mimoChan.indGOP = 1;
    mimoChan.indCoherent = 1;
end
