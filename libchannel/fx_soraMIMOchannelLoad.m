function [H, N, Tc] = fx_soraMIMOchannelLoad(localCSIpath, SORAscen, maxNr)
    
    global mimoChanfileCSI;
    global mimoChanfileInd;
    global mimoChanpktInd;
    global mimoChanFreqInd;
    
    SORAscen = SORAscen(1:5);
    
    Tc = 152 / maxNr;
    flagFoundData = false;
    while ~flagFoundData
        if isempty(mimoChanfileCSI)
            matfile = sprintf('%s/%s_%d.mat', localCSIpath, SORAscen, mimoChanfileInd);
            A = load(matfile);
            if isempty(A.fileCSI)
                mimoChanfileInd = mimoChanfileInd + 1;
            else
                mimoChanfileCSI = A.fileCSI;
                mimoChanpktInd = 1;
                mimoChanFreqInd = 1;
            end
        else
            if numel(mimoChanfileCSI{mimoChanpktInd}.H) < maxNr
                %check the number of working radio packet by packet
                mimoChanpktInd = mimoChanpktInd + 1;
                if mimoChanfileInd == 32 && mimoChanpktInd > numel(mimoChanfileCSI)
                    error('radio fail');
                end
            else
                % load data frequency by frequency
                H = [];
                N = [];
                for nr = 1:numel(mimoChanfileCSI{mimoChanpktInd}.H)
                    H = [H; mimoChanfileCSI{mimoChanpktInd}.H{nr}(mimoChanFreqInd,:)];
                    N = [N; mimoChanfileCSI{mimoChanpktInd}.N{nr}(mimoChanFreqInd,:)];
                end
                mimoChanpktInd = mimoChanpktInd + floor((mimoChanFreqInd + 1) / 49);
                mimoChanFreqInd = mod(mimoChanFreqInd, 48) + 1;
                flagFoundData = true;
            end
            if mimoChanpktInd > numel(mimoChanfileCSI)
                mimoChanfileInd = mimoChanfileInd + 1;
                mimoChanfileCSI = [];
            end
        end
        if mimoChanfileInd > 32
            error('touch the end');
        end
    end
end