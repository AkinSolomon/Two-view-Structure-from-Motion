function error = reproj2(p, mp, npics, numpts)
 
    M = zeros(3,4,npics);
    pointStart = 3*npics+1;
    P = mp(:,pointStart:end);
    errorC=[];
    error=zeros(1,2);
    for i=1:npics
        split = 3*(i-1)+1;
        M(:,:,i) = mp(:,split:2+split)';
        i;
        for j=1:numpts
            error(1) = x(p,i,j) - (M(1,:,i)*P(:,j))/(M(3,:,i)*P(:,j));
            error(2) = y(p,i,j) - (M(2,:,i)*P(:,j))/(M(3,:,i)*P(:,j));
            errorC=[errorC,error'];
        end
    end
    error=errorC;
    %error = sum(errorTemp)/numPics/numPoints;
end

function xV = x(p,i,j)
    xV=p(3*(i-1)+1,j);
end

function yV = y(p, i,j)
    yV=p(3*(i-1)+2,j);
end