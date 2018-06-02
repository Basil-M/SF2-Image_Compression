function [vlc, bits, huffval] = dwt_enc(X,N_LEVELS, M,q0,rise, N_sup,N_LBT, opthuff, dcbits)
    
% DWT_ENC Encodes an image to a (simplified) JPEG bit stream
%
%  [vlc bits huffval] = jpegenc(X, qstep, N, M, opthuff, dcbits) Encodes the
%  image in X to generate the variable length bit stream in vlc.
%
%  Yq is the encoded output of either a DWT/LBT/DCT
%  N is the width of DCT block OR number of levels of DWT
%  subim_width is the width of each block to be coded (defaults to N). Must be an
%  integer multiple of N - if it is larger, individual blocks are
%  regrouped.
%  if original data X is provided, it will optimise the code based on it
%  dcbits determines how many bits are used to encode the DC coefficients
%  of the DCT (defaults to 8)
%
%  vlc is the variable length output code, where vlc(:,1) are the codes, and
%  vlc(:,2) the number of corresponding valid bits, so that sum(vlc(:,2))
%  gives the total number of bits in the image
%  bits and huffval are optional outputs which return the Huffman encoding
%  used in compression

% This is global to avoid too much copying when updated by huffenc
global huffhist  % Histogram of usage of Huffman codewords.

% Presume some default values if they have not been provided
if ((nargout~=1) && (nargout~=3)) error('Must have one or three output arguments'); end
if ~exist('dcbits','var')
    dcbits = 16;
end
if ~exist('opthuff')
    opthuff = true;
elseif opthuff == -1;
    opthuff = true;
end

if ~exist('M','var')
    M = 2^N_LEVELS;
end

if ~exist('rise','var')
    rise = 1;
end

if ~exist('N_sup','var')
    N_sup = -1;
end

if ~exist('N_LBT','var')
    N_LBT = -1;
end
% make it zero mean
%if sum(X(:) > 128) ~= 0; X  = X - 128; end
%if ((opthuff==true) && (nargout==1)) error('Must output bits and huffval if optimising huffman tables'); end
N = 2^N_LEVELS;
% DWT on input image X.
%fprintf(1, 'Forward %i level DWT\n', N_LEVELS);
Y = nlevdwt(X, N_LEVELS);

if N_sup > 0
    Y = dwt_suppress(Y, N_LEVELS, N_sup);
end
q_rats = dwt_q_ratios(size(Y),N_LEVELS);

if N_LBT > 0
    q_rats(1,N_LEVELS+1) = q_rats(1,N_LEVELS+1);
    m_lbt = length(Y)/2^N_LEVELS;
    q_lbt = q0*q_rats(1, N_LEVELS+1);
    Y(1:m_lbt, 1:m_lbt) = lbt_enc(Y(1:m_lbt, 1:m_lbt), 4, q_lbt, sqrt(2), rise);
end

% Quantise to integers.
%fprintf(1, 'Quantising DWT', q0); 
Yq=dwtquant1(Y,N_LEVELS, q0*q_rats,rise);

% reshuffle to look like DCT
%fprintf(1, 'Regrouping %i level DWT to look like %i x %i DCT', N_LEVELS, N, N);
Yq = dwtgroup(Yq, N_LEVELS);
% Generate zig-zag scan of AC coefs.
scan = diagscan(M);

% On the first pass use default huffman tables.
%disp('Generating huffcode and ehuf using default tables')
[dbits, dhuffval] = huffdflt(1);  % Default tables.
[huffcode, ehuf] = huffgen(dbits, dhuffval);

% Generate run/ampl values and code them into vlc(:,1:2).
% Also generate a histogram of code symbols.
%disp('Coding rows')
sy=size(Yq);
t = 1:M;
huffhist = zeros(16*16,1);
vlc = [];
for r=0:M:(sy(1)-M),
  vlc1 = [];
  for c=0:M:(sy(2)-M),
    yq = Yq(r+t,c+t);
    % Possibly regroup 
    if (M > N) yq = regroup(yq, N); end
    % Encode DC coefficient first
    yq(1) = yq(1) + 2^(dcbits-1);
    if ((yq(1)<1) | (yq(1)>(2^dcbits-1)))
      error('DC coefficients too large for desired number of bits');
    end
    dccoef = [yq(1)  dcbits]; 
    % Encode the other AC coefficients in scan order
    ra1 = runampl(yq(scan));
    vlc1 = [vlc1; dccoef; huffenc(ra1, ehuf)]; % huffenc() also updates huffhist.
  end
  vlc = [vlc; vlc1];
end

% Return here if the default tables are sufficient, otherwise repeat the
% encoding process using the custom designed huffman tables.
if (opthuff==false) 
  if (nargout>1)
    bits = dbits;
    huffval = dhuffval;
  end
  %fprintf(1,'Bits for coded image = %d\n', sum(vlc(:,2)));
  return;
end

% Design custom huffman tables.
%disp('Generating huffcode and ehuf using custom tables')
[dbits, dhuffval] = huffdes(huffhist);
[huffcode, ehuf] = huffgen(dbits, dhuffval);

% Generate run/ampl values and code them into vlc(:,1:2).
% Also generate a histogram of code symbols.
%disp('Coding rows (second pass)')
t = 1:M;
huffhist = zeros(16*16,1);
vlc = [];
for r=0:M:(sy(1)-M),
  vlc1 = [];
  for c=0:M:(sy(2)-M),
    yq = Yq(r+t,c+t);
    % Possibly regroup 
    if (M > N) yq = regroup(yq, N); end
    % Encode DC coefficient first
    yq(1) = yq(1) + 2^(dcbits-1);
    dccoef = [yq(1)  dcbits]; 
    % Encode the other AC coefficients in scan order
    ra1 = runampl(yq(scan));
    vlc1 = [vlc1; dccoef; huffenc(ra1, ehuf)]; % huffenc() also updates huffhist.
  end
  vlc = [vlc; vlc1];
end
%fprintf(1,'Bits for coded image = %d\n', sum(vlc(:,2)))
%fprintf(1,'Bits for huffman table = %d\n', (16+max(size(dhuffval)))*8)

if (nargout>1)
  bits = dbits;
  huffval = dhuffval';
end
return
