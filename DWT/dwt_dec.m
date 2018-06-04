function Z = dwt_dec(vlc, N_LEVELS,M, q0,rise,N_LBT, bits, huffval, dcbits, W, H)

% JPEGDEC Decodes a (simplified) JPEG bit stream to an image
%
%  Z = jpegdec(vlc, qstep, N, M, bits huffval, dcbits, W, H) Decodes the
%  variable length bit stream in vlc to an image in Z.
%
%  vlc is the variable length output code from jpegenc
%  qstep is the quantisation step to use in decoding
%  N is the width of the DCT block (defaults to 8)
%  M is the width of each block to be coded (defaults to N). Must be an
%  integer multiple of N - if it is larger, individual blocks are
%  regrouped.
%  if bits and huffval are supplied, these will be used in Huffman decoding
%  of the data, otherwise default tables are used
%  dcbits determines how many bits are used to decode the DC coefficients
%  of the DCT (defaults to 8)
%  W and H determine the size of the image (defaults to 256 x 256)
%
%  Z is the output greyscale image

% Presume some default values if they have not been provided
dwt_scan = (M < 1);
opthuff = exist('bits','var')&exist('huffval','var');
if ~exist('H','var') || ~exist('W','var')
    H = 256; W = 256;
end
if ~exist('dcbits','var')
    dcbits = 16;
end

%DWT is reshuffled to look like DCT of size 2^N
N_DCT = 2^N_LEVELS;

if ~exist('M', 'var')
    M = N_DCT;
elseif M == -1
    M = N_DCT;
end



% Set up standard scan sequence
if dwt_scan
    scan = dwtscan(W, N_LEVELS);
else
    if (mod(M, N_DCT)~=0) error('Encoding width must be an integer multiple of 2^N_LEVELS'); end
    scan = diagscan(M);
end

if (opthuff)
  %disp('Generating huffcode and ehuf using custom tables')
else
  %disp('Generating huffcode and ehuf using default tables')
  [bits, huffval] = huffdflt(1);
end
% Define starting addresses of each new code length in huffcode.
huffstart=cumsum([1; bits(1:15)]);
% Set up huffman coding arrays.
[huffcode, ehuf] = huffgen(bits, huffval);

% Define array of powers of 2 from 1 to 2^16.
k=[1; cumprod(2*ones(16,1))];

% For each block in the image:

% Decode the dc coef (a fixed-length word)
% Look for any 15/0 code words.
% Choose alternate code words to be decoded (excluding 15/0 ones).
% and mark these with vector t until the next 0/0 EOB code is found.
% Decode all the t huffman codes, and the t+1 amplitude codes.

eob = ehuf(1,:);
run16 = ehuf(15*16+1,:);
i = 1;
Zq = zeros(H, W);
t=1:M;
if ~dwt_scan
%disp('Decoding rows')
    for r=0:M:(H-M),
      for c=0:M:(W-M),
        yq = zeros(M,M);

    % Decode DC coef - assume no of bits is correctly given in vlc table.
        cf = 1;
        if (vlc(i,2)~=dcbits) error('The bits for the DC coefficient does not agree with vlc table'); end
        yq(cf) = vlc(i,1) - 2^(dcbits-1);
        i = i + 1;


    % Loop for each non-zero AC coef.
        while any(vlc(i,:) ~= eob),
          run = 0;

    % Decode any runs of 16 zeros first.
          while all(vlc(i,:) == run16), run = run + 16; i = i + 1; end

    % Decode run and size (in bits) of AC coef.
          start = huffstart(vlc(i,2));
          res = huffval(start + vlc(i,1) - huffcode(start));
          run = run + fix(res/16);
          cf = cf + run + 1;  % Pointer to new coef.
          si = rem(res,16);
          i = i + 1;

    % Decode amplitude of AC coef.
          if vlc(i,2) ~= si,
            error('Problem with decoding .. you might be using the wrong bits and huffval tables');
            return
          end
          ampl = vlc(i,1);

    % Adjust ampl for negative coef (i.e. MSB = 0).
          thr = k(si);
          yq(scan(cf-1)) = ampl - (ampl < thr) * (2 * thr - 1);

          i = i + 1;      
        end

    % End-of-block detected, save block.
        i = i + 1;

        % Possibly regroup yq
        if (M > N_DCT) yq = regroup(yq, M/N_DCT); end
        Zq(r+t,c+t) = yq;
      end
    end
else
    %encode DC coefficients
    m_lp = W/2^N_LEVELS;    %size of the DC portion 
    i = 1;
    cf = 1;
    for r = 1: m_lp
        for c = 1:m_lp
            if (vlc(i,2)~=dcbits) error('The bits for the DC coefficient does not agree with vlc table'); end
            Zq(r,c) = vlc(i, 1) - 2^(dcbits - 1);
            i = i + 1;
        end
    end
    % Loop for each non-zero AC coef.
        while(i <= length(vlc))
          run = 0;

    % Decode any runs of 16 zeros first.
          while all(vlc(i,:) == run16), run = run + 16; i = i + 1; end

    % Decode run and size (in bits) of AC coef.
          if(i <= length(vlc))
              start = huffstart(vlc(i,2));
              res = huffval(start + vlc(i,1) - huffcode(start));
              run = run + fix(res/16);
              cf = cf + run + 1;  % Pointer to new coef.
              si = rem(res,16);
              i = i + 1;
          end

    % Decode amplitude of AC coef.
          if(i <= length(vlc))
              if vlc(i,2) ~= si,
                error('Problem with decoding .. you might be using the wrong bits and huffval tables');
              end
              ampl = vlc(i,1);

        % Adjust ampl for negative coef (i.e. MSB = 0).
              thr = k(si);
              Zq(scan(cf-1)) = ampl - (ampl < thr) * (2 * thr - 1);
          end
          i = i + 1;      
        end
end

if isscalar(q0)
    q = q0*dwt_q_ratios(size(Zq), N_LEVELS);
else
    q = q0;
end

if N_LBT > -1
    q(1,N_LEVELS+1) = 0;%q(1,N_LEVELS+1)*1;
end

if ~exist('rise')
    rise = 1;
end

if(~dwt_scan)
    fprintf(1, 'Unregrouping DWT\n');
    Zq = dwtgroup(Zq, -1*N_LEVELS);
end
fprintf('Inverse quantising the DWT');
Zi=dwtquant2(Zq,N_LEVELS,q, rise);  


fprintf(1, 'Inverse %i level DWT\n', N_LEVELS);


if N_LBT > 0
    disp('Inverse LBT on lowpass component');
    m_lbt = length(Zq)/2^N_LEVELS;
    q_lbt = q(1, N_LEVELS+1);
    N_LBT = min(m_lbt, 4);
    Zi(1:m_lbt, 1:m_lbt) = lbt_dec(Zi(1:m_lbt, 1:m_lbt), N_LBT, sqrt(2));
end


Z = nlevidwt(Zi, N_LEVELS);

%if(sum(Z(:) > 128)) == 0
Z = Z + 128;
%end
return