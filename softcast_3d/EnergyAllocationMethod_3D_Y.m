function [ G_Y ] = EnergyAllocationMethod_3D_Y(TxEnergy, VarVid_Y)
%
G_Y = zeros(size(VarVid_Y));
G_Y(VarVid_Y ~= 0) = (VarVid_Y(VarVid_Y ~= 0)) .^ (-0.25);

EnergyTotal = numel(G_Y) * TxEnergy;
EnergyPreScale = sum(sum(sum((G_Y .^ 2) .* VarVid_Y)));
ScaleFactor = (EnergyPreScale / EnergyTotal) .^ -0.5;

G_Y = G_Y .* ScaleFactor;
end