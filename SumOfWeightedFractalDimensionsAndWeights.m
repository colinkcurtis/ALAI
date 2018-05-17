function [y,z] = SumOfWeightedFractalDimensionsAndWeights ...
	(i, AllRoughnessParameters)
    
    CurrentDf = AllRoughnessParameters{1,i}{2,1}(1,1);
    
    CurrentDfError = AllRoughnessParameters{1,i}{3,1}(1,1);
    
    CurrentDfWeight = 1 / (CurrentDfError)^2; 
    
    y = CurrentDf*CurrentDfWeight;
    
    z = CurrentDfWeight;
    