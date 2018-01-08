function [P, M] = triangulate(points, K, numPoints,k,l)
    
    % P is 3D points
    % M is camera matrix
    % k is the first image
    % l is the second image

    p1 = points(3*k-2:3*k,:);
    p2 = points(3*l-2:3*l,:);
    
    %map points to back to camera plane
    p1=inv(K)*p1;
    p2=inv(K)*p2;

    F = myfund(p1, p2);

    [U, ~, V] = svd(F);
    sk1 = [0 1 0 ; -1 0 0 ; 0 0 1];
    sk2 = [0 1 0 ; -1 0 0 ; 0 0 0];
    
    Rp = U*sk1*-V';
    Ts = U*sk2'*U';
    Tp = [Ts(3,2);Ts(1,3);Ts(2,1)];
    M2 = [Rp, Tp];
    M1 = eye(3,4);

    for j = 1:numPoints
        A = [p1(1,j)*M1(3,:)- M1(1,:); p1(2,j)*M1(3,:) - M1(2,:);
              p2(1,j)*M2(3,:) - M2(1,:); p2(2,j)*M2(3,:) - M2(2,:)];
        [~, ~, V] = svd(A'*A);
        P(:,j) = V(:,4) / V(4,4);
    end
    
    M=M2;
end
