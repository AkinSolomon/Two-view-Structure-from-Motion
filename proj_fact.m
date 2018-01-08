
%Returns the SVD factorization of a matrix of observations, W
%
function [M P]=proj_fact(W)

m = size(W,1)/3;
   n = size(W,2);
   
   D = ones(m,n);

   for it = 1:2
       for p = 1:n
           lp = D(:,p);
           D(:,p) = lp/norm(lp);
       end
       for i = 1:m
           li = D(i,:);
           D(i,:) = li/norm(li);
      end
   end
   
   for i = 1:m
      for p = 1:n
          W(3*i-2:3*i,p) = D(i,p)*W(3*i-2:3*i,p);
      end
   end


[U S V]=svd(W);

[U,S,V] = svd(W);
S = S/S(1,1);
M=U(:,1:4);
P=S(1:4,1:4)*V(:,1:4)';

end