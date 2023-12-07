function [X,IX]=PROX(x,X,U,a)
%U=roundn(U,-4);
S=U(:);
[~,IX] = sort(S(:),'descend');
if(a>length(IX))
    a=length(IX);
end
x(IX(1:a))=U(IX(1:a));
% x=reshape(x,size(X,1),size(X,2));
X=x;
X(X<0)=0;
end

