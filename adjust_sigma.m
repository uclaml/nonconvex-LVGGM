function [hsigma] = adjust_sigma(hsigma)
d = size(hsigma,1);
eigVal = eig(hsigma);
minV = min(eigVal);
maxV = max(eigVal);
alpha = 1/max(abs(maxV),abs(minV));
hsigma = (1-alpha*abs(minV))*hsigma + abs(minV)*eye(d);
end

