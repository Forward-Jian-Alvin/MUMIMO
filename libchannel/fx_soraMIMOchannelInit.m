function fx_soraMIMOchannelInit(SORApath, SORAscen, SORAinfopath, localCSIpath)

    global mimoChanfileInd;
    global mimoChanpktInd;
    global mimoChanFreqInd;

    datamap = {'//msraim-clu16/SORA_MIMO_CSI', ...
               '//msraim-clu15/SORA_MIMO_CSI', ...
               '//msraim-clu14/SORA_MIMO_CSI', ...
               '//msraim-clu13/SORA_MIMO_CSI', ...
               '//msraim-clu12/SORA_MIMO_CSI', ...
               '//msraim-clu09/SORA_MIMO_CSI', ...
               '//msraim-clu08/SORA_MIMO_CSI', ...
               '//msraim-clu07/SORA_MIMO_CSI'};
    if ~isempty(regexp(SORAscen, 'D\d_\d\d-\d', 'once'))
        indState = sscanf(SORAscen(6:end), '-%d');
    else
        indState = 0;
    end
	SORAscen = SORAscen(1:5);
    inddata = sscanf(SORAscen, 'D%d');
    
    if exist(localCSIpath, 'dir') == 0
        mkdir(localCSIpath);
    end
    for ii = 1:32
        matfile = sprintf('%s/%s_%d.mat', localCSIpath, SORAscen, ii);
        
        if exist(matfile, 'file') == 0
            srcfile = sprintf('%s/%s_%d.mat', datamap{inddata}, SORAscen, ii);
            if exist(srcfile, 'file') ~= 0
                copyfile(srcfile, matfile);
            else
                channeldir = sprintf('%s/%s/%d/', SORApath, SORAscen, ii);
                infofile = sprintf('%s/%s/traceinfo_%s_%d.txt', SORAinfopath, SORAscen, SORAscen, ii);
                info = load(infofile);

                fileCSI = cell(size(info,1)-1,1);
                for jj = 1:size(info,1)-1 %the last line is the hist info of working radio
                    pktCSI = struct('H', cell(1,1), 'N', cell(1,1));
                    ind = 1;
                    for kk = 0:3
                        if info(jj,3+kk) ~= -100
                            Hfile = sprintf('%s/%d_rad%d_P%d_S%d_H.txt', channeldir, ii, kk, info(jj,3+kk), info(jj,7+kk));
                            Nfile = sprintf('%s/%d_rad%d_P%d_S%d_N.txt', channeldir, ii, kk, info(jj,3+kk), info(jj,7+kk));
                            chanH = load(Hfile);
                            chanN = load(Nfile);

                            pktCSI.H{ind} = chanH(:,1:2:end) + 1i*chanH(:,2:2:end);
                            pktCSI.N{ind} = chanN(:,1:2:end) + 1i*chanN(:,2:2:end);
                            ind = ind + 1;
                        end
                    end
                    fileCSI{jj} = pktCSI;
                end
                save(matfile, 'fileCSI');
            end
        end
    end
    
    if indState > 0
        StateInfoFile = sprintf('SNR-%s.txt', SORAscen);
        localStateInfoFile = sprintf('%s/%s', localCSIpath, StateInfoFile);
        if exist(localStateInfoFile, 'file') == 0
            copyfile([datamap{inddata} '/' StateInfoFile], localStateInfoFile);
        end
        A = importdata(localStateInfoFile);
        mimoChanfileInd = A(indState, 3);
        mimoChanpktInd  = A(indState, 4);
        mimoChanFreqInd = A(indState, 5);
    else
        fx_soraMIMOchannelReset();
    end
    fx_soraMIMOchannelSave(1);
end