function out = compute_norms(streamer_out)
errors = streamer_out(2);
errors = errors{1};
pnerrors = errors(1,:,1,1,end);
l1 = sum(abs(pnerrors))/size(pnerrors,2);
l2 = sqrt(sum(pnerrors.^2)/size(pnerrors,2));
linf = max(abs(pnerrors));
out = [l1, l2, linf];
end