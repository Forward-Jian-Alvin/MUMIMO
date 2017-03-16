function H = fx_MIMOchannelLoad(Nr)

global mimoChan;

% debug begin
% fprintf(1, '%d, %d\n',mimoChan.indFrm, mimoChan.indCoherent);
% debug end

H = mimoChan.Hrayleigh(1:Nr, :, mimoChan.indGOP, mimoChan.indCoherent);

mimoChan.indCoherent = mimoChan.indCoherent + 1;
