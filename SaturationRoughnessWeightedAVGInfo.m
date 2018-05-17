function y = SaturationRoughnessWeightedAVGInfo(SumOfWeightedSatRuffs, ...
	SumOfSatRuffWeights, SatRuffWeightedAVGError)
    
    SatRuffWeightedAVG = SumOfWeightedSatRuffs / SumOfSatRuffWeights;
    
	SatRuffWeightedAVGError = ceil(SatRuffWeightedAVGError/10^floor(log10 ...
		(SatRuffWeightedAVGError)))...
			*10^floor(log10(SatRuffWeightedAVGError)); 
			% salt to suit...   
    
	SatRuffErrorInfoAVG = num2str(SatRuffWeightedAVGError);
    
	SatRuffErrorInfoAVG = cat(2,'(',SatRuffErrorInfoAVG,')');
    
	SatRuffWeightedAVG = num2str(SatRuffWeightedAVG);
    
	SatRuffWeightedAVG = SatRuffWeightedAVG(1:end-2);
    
	y = 'Saturation Roughness: ';
    
	y = cat(2,y, SatRuffWeightedAVG,' ',SatRuffErrorInfoAVG,' [A]');