% works after the run, to get a histogram of Ho values
% hard coded to work on i_density_rad15_step2.mat

load i_density_rad15_step2.mat

v = 4*pi*radius^3 / 3;
tmp = maxd.*v;
figure;hist(tmp,unique(tmp))
p = prctile(tmp,95);
tmp2 = tmp(tmp >= p);
hold on; h = hist(tmp2,unique(tmp2));
bar(unique(tmp2),h,'r')

for i = p+1:-1:1
    pp(i) = sum(tmp > i) ./ length(tmp);
    xx(i) = i;
end
figure('Color','w'),plot(xx,pp,'rs-','LineWidth',2)
xlabel('Number of peaks within radius'),ylabel('Corrected p')

disp(['Threshold is at least ' num2str(p) ' points within ' num2str(radius) ' mm'])
