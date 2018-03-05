function Y = HardThresh(X,s)
d = size(X,1);
indx = round((s-d)/2);

X_tri = triu(X);
X_diag = diag(diag(X));
X_sorted = X_tri - X_diag;
X_sorted = sort(abs(X_sorted(:)),'descend');

Y_tri = wthresh(X_tri-X_diag,'h',X_sorted(indx));
Y = Y_tri + Y_tri' + X_diag;
end