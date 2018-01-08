function [M,Points] = splitMP(MP, numPics)

    M = zeros(3,4,numPics);
    pointStart = 3*numPics+1;
    Points = MP(:,pointStart:end);
    errorC=[];
    error=zeros(1,2);
    
    for i=1:numPics
        split = 3*(i-1)+1;
        M(:,:,i) = MP(:,split:2+split)';
    end
    
    error=errorC;
    %error = sum(errorTemp)/numPics/numPoints;
    
end