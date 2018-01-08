function [ F ] = myfund( p1, p2)
%motivated by Taehee Lee's vision toolbox
n = size(p1, 2);

[p1,H1]=normalise2dpts(p1);
[p2,H2]=normalise2dpts(p2);

%solving x'Fx=0 using all the pairs of points

X = zeros(size(p1,2), 9 );
for i = 1:n
    X(i,:) = kron(p1(:,i), p2(:,i));
end

[~,~, Vx] = svd(X);

%smallest singular value of V
Fs = Vx(:,9);

% unstack F
F = reshape(Fs, 3, 3);

[U D V] = svd(F);
s = (D(1,1)+D(2,2))/2;
D = diag([ s s 0 ]);

F = U * D * V';

% apply normalization of x1 and x2
F = H2' * F * H1;

end
