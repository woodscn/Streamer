function out = fit_errors(dxes,norms)
l1s = zeros(1,size(norms,2));
l2s = zeros(1,size(norms,2));
linfs = zeros(1,size(norms,2));
for ind = 1:size(norms,2)
    l1s(ind) = norms(ind,1);
    l2s(ind) = norms(ind,2);
    linfs(ind) = norms(ind,3);
end
p2 = polyfit(log(dxes),log(l2s),1);
coefficient = exp(p2(2));
exponent = p2(1);
out = [coefficient, exponent];
end