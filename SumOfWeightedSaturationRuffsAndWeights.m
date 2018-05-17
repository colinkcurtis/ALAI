function [y,z] = SumOfWeightedSaturationRuffsAndWeights ...
	(i, AllRoughnessParameters)
    
    CurrentSatRuff = AllRoughnessParameters{1,i}{2,1}(1,2);
    
    CurrentSatRuffError = AllRoughnessParameters{1,i}{3,1}(1,2);
    
    CurrentSatRuffWeight = 1 / (CurrentSatRuffError)^2;
    
    y = CurrentSatRuff*CurrentSatRuffWeight;
    
    z = CurrentSatRuffWeight;