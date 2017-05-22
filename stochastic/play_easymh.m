function [smpl] = mh(x0, nsteps, pdf)
% using normal distr as jumping fn

x = x0(:)';
smpl = zeros(nsteps, size(x,2));
for i = 1:nsteps
    candidate = x + randn; % or randn(size(x))??
    acceptance = min(1, pdf(candidate)/pdf(x));
    if rand < acceptance
        x = candidate;
    end
    smpl(i,:) = x;
end



% % trying my mh.m function
% x0 = 1;
% nsteps = 15000;
% pdf = @(x) normpdf(x);
% smpl = mh(x0, nsteps, pdf);
% figure;
% h = histfit(smpl,100);
% h(1).FaceColor = [.8 .8 1];
