%CDS Matlab script 
%2014 Georg Ramer (georg@ramer.at)
%implementing the description of CDS given in
%
%I. Noda, “Two-dimensional codistribution spectroscopy to determine the 
%   sequential order of distributed presence of species,” 
%   J. Mol. Struct., pp. 1–10, Jan. 2014.
%





function [ sync, async ] = CDS_cross( data1, data2 )
%CDS calculates the codistribution spectra for the given data
% spectra are entered one spectrum per column

    data = vertcat(data1, data2);
    
    rows_data1 = (1:size(data,1)) <= size(data1,1);
    
    meancentered = zeros(size(data));
    rowmean = mean(data, 2);
    for col = 1:size(data,2)
        meancentered(:, col) = data(:, col) - rowmean;
    end
    
    
    autocorr = sum(meancentered .* meancentered,2);
    T = autocorr * autocorr.';
    m = size(data,2);
    ks = 1: (m);
    summands = sum(row_multiply( meancentered, ks) ,2) ./ rowmean;
    summands = repmat(summands,  1, length(summands));
    async = (summands - summands') .* T / (m * (m-1));
    
    sync = sqrt(T.^2 - async.^2);  
    
    
    sync = sync(~rows_data1, rows_data1);
    async = async (~rows_data1, rows_data1);
end

function [res] = row_multiply(data, row)
    res = zeros(size(data));
    for r = 1:size(data, 1)
        res = data(r,:) .* row;
    end
end

