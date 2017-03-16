function fx_soraMIMOchannelSet(indGOP, cfg)
    global State;
    global mimoChanfileInd;
    global mimoChanpktInd;
    global mimoChanFreqInd;
    global mimoChanfileCSI;
    
    if ~isempty(State(indGOP))
        mimoChanfileInd = State{indGOP}(1);
        mimoChanpktInd  = State{indGOP}(2);
        mimoChanFreqInd = State{indGOP}(3);

        SORAscen =  cfg.SORAscen(1:5);
        matfile = sprintf('%s/%s_%d.mat', cfg.localCSIpath, SORAscen, mimoChanfileInd);
        A = load(matfile);
        mimoChanfileCSI = A.fileCSI;
    end
end