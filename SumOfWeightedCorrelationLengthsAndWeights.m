function [y,z] = SumOfWeightedCorrelationLengthsAndWeights ...
	(i, AllRoughnessParameters)
        
    CurrentCorrLength = AllRoughnessParameters{1,i}{2,1}(1,3);
    
    CurrentCorrLengthError = AllRoughnessParameters{1,i}{3,1}(1,3);
    
    CurrentCorrLengthWeight = 1 / (CurrentCorrLengthError)^2;
    
    y = CurrentCorrLength*CurrentCorrLengthWeight;
    
    z = CurrentCorrLengthWeight;