% numerical integral of the overlapping area of two normal distributions:
% s1,s2...sigma of the normal distributions 1 and 2
% mu1,mu2...center of the normal distributions 1 and 2
% xstart,xend,xinterval...defines start, end and interval width
% example: [overlap] = calc_overlap_twonormal(2,2,0,1,-10,10,0.01)
function [overlap2] = calc_overlap_twonormal(s1,s2,mu1,mu2,xstart,xend,xinterval)
%clf
x_range=xstart:xinterval:xend;
%plot(x_range,[normpdf(x_range,mu1,s1)' normpdf(x_range,mu2,s2)']);
%hold on
%area(x_range,min([normpdf(x_range,mu1,s1)' normpdf(x_range,mu2,s2)']'));
overlap=cumtrapz(x_range,min([normpdf(x_range,mu1,s1)' normpdf(x_range,mu2,s2)']'));
overlap2 = overlap(end);
legend([num2str(overlap2)]);