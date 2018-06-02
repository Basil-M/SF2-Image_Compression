function [scan] = dwtscan(S,N)

% DWTSCAN Generates scanning pattern for n-levelDWT of S x S image
%
%  [scan] = DIAGSCAN(N) Produces a diagonal scanning index for
%  an NxN matrix
%
%  The first entry in the matrix is assumed to be the DC coefficient
%  and is therefore not included in the scan
v = @(row,col) S*(row-1) + col; %convert (r,c) into column vector

s2 = S/2;
%start from top left corner of vertical component
scan = [];
%vertical scan
for i = 1:N
    for c = (1+s2):S
        for r = 1:s2
            scan = [scan v(r,c)];
        end
    end
    %horizontal scan
    for r = (1+s2):S
        for c = 1:s2
            scan = [scan v(r,c)];
        end
    end
    %diagonal scan on diagonal details
    scan_diag = [1 diagscan(s2)];
    rs = ceil(scan_diag/s2);
    scan_diag = s2*S + rs*s2 + scan_diag;
    scan = [scan scan_diag];
    
    S = S/2;
    s2 = s2/2;
end
end

