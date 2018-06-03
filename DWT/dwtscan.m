function [scan] = dwtscan(S,N)

% DWTSCAN Generates scanning pattern for n-levelDWT of S x S image
%
%  [scan] = DIAGSCAN(N) Produces a diagonal scanning index for
%  an NxN matrix
%
%  The first entry in the matrix is assumed to be the DC coefficient
%  and is therefore not included in the scan
v = @(row,col) S*(row-1) + col; %convert (r,c) into column vector
S_c = S;
s2 = S_c/2;
%start from top left corner of vertical component
scan = [];
%vertical scan
for i = 1:N
    for c = (1+s2):2:S_c
        for r = 1:s2
            scan = [scan v(r,c)];
        end
        if c~=S_c
            %and back up
            for r = s2:-1:1
                scan = [scan v(r,c+1)];
            end
        end
    end
    %horizontal scan
    for r = (1+s2):2:S_c
        for c = 1:s2
            scan = [scan v(r,c)];
        end
        if r~=S_c
            %and back to the left
            for c = s2:-1:1
                scan = [scan v(r+1,c)];
            end
        end
    end
    %diagonal scan on diagonal details
    scan_diag = [1 diagscan(s2)];
    s_convert = s2*S + s2;      %coordinates of start point of this sub square
                                %in global square
    
    rs = ceil(scan_diag/s2)-1;   %number of rows from that start point of each point
    
    scan_diag = rs*(S - s2) + s_convert + scan_diag;
    scan = [scan scan_diag];
    
    S_c = S_c/2;
    s2 = s2/2;
end
end

