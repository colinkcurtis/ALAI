function [y,z] = ImageSizeProperties(Pixels, AFMMetaData)
    
    tempStringforPhysicalImageSize= AFMMetaData{1}{1}{1};
    clippedString = strrep(tempStringforPhysicalImageSize,'ScanSize: ','');
    clippedString1 = 10e9*str2num(clippedString); 
    clippedString = num2str(clippedString1);
    y = cat(2, 'Scan Size: ', ' ', clippedString, ' [A]'); 
    
    ImageResolution = clippedString1/Pixels(1,1);
    ImageResolution = ceil(ImageResolution/10^floor ...
		(log10(ImageResolution)))*10^floor(log10(ImageResolution));
    ImageResolution = num2str(ImageResolution);
    z = cat(2, 'Image Resolution: ', ImageResolution, ' [A/pixel]');