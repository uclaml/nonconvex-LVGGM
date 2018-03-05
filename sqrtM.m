function B = sqrtM(A)
    [V, D,U] = svd(A);
    diagD = diag(D);
    index = find(diagD>=0);
    sqrtD = diag(sqrt(diagD(index)));
    B = V(:,index)*sqrtD*U(:,index)';
end