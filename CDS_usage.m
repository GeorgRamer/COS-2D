% How to use the CDS functions

load('reference_data.mat')

figure()
plot(wn, data);

figure()
surf(1:size(data,2), wn, data, 'EdgeColor', 'none')


%% one wavelength axis

[sync,async] = CDS(data);

COS_plot(wn, wn, sync, 10)
title('synchronous plot')



COS_plot(wn,wn,async,10)
title('asynchronous plot')



