function fx_soraMIMOchannelSave(indGOP)
    
    global mimoChanfileInd;
    global mimoChanpktInd;
    global mimoChanFreqInd;
    
    global State;
    State{indGOP} = [mimoChanfileInd mimoChanpktInd mimoChanFreqInd];
end