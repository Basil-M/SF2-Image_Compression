function Z = jpegdeclbtnlev(vlc, qstep, N, M, rise1, s, bits, huffval, dcbits, ratio, ratio2, W, H)
% performing a multi level lbt

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
%  s is the step size for the LBT. This defaults to sqrt(2)
%  if bits and huffval are supplied, these will be used in Huffman decoding
%  of the data, otherwise default tables are used
%  dcbits determines how many bits are used to decode the DC coefficients
%  of the DCT (defaults to 8)
%  ratio is the ratio between the 1st level of the lbt and the second level
%  W and H determine the size of the image (defaults to 256 x 256)
%
%  Z is the output greyscale image

% Presume some default values if they have not been provided
error(nargchk(2, 13, nargin, 'struct'));
opthuff = true;
if (nargin<13)
  H = 256;
  W = 256;
    if(nargin<11)
        ratio2 = 0;
      if(nargin<10)
          ratio=1;
          if (nargin<9)
            dcbits = 16; % changed
            if (nargin<8)
              opthuff = false;
              if(nargin<6)
                  s = sqrt(2);
                  if(nargin<5)
                      rise1 = 1;
                      if (nargin<4)
                        if (nargin<3)
                          N = 8;
                          M = 8;
                        else
                          M = N;
                        end
                      else 
                        if (mod(M, N)~=0) error('M must be an integer multiple of N'); end
                      end
                  end
              end
            end
          end
     end
  end
end

% Set up standard scan sequence
scan = diagscan(M);

if (opthuff)
  disp('Generating huffcode and ehuf using custom tables')
else
  disp('Generating huffcode and ehuf using default tables')
  [bits huffval] = huffdflt(1);
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

disp('Decoding rows')
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
    if (M > N) yq = regroup(yq, M/N); end
    Zq(r+t,c+t) = yq;
  end
end

fprintf(1, 'Inverse quantising to step size of %i\n', qstep);
%Work out how to remove the quantisation when I'm quantising each bit
%separately
Zi=quant2(Zq,qstep,rise1*qstep);
Zi(1:N:W, 1:N:H)=quant2(Zq(1:N:W, 1:N:H),ratio*qstep, ratio*qstep*rise1);
Zi(1:N^2:W, 1:N^2:H) = quant2(Zq(1:N^2:W,1:N^2:H),ratio2*ratio*qstep,ratio2*ratio*qstep*rise1);

% add in if statement for this
    % decode the encoded 3rd level dc coefficients
    %Zj = Zi(1:N:W, 1:N:H);
    %Zj(1:N:W/N,1:N:H/N) = lbt_dec(Zj(1:N:W/N,1:N:H/N),N,s)*N;
    %Zi(1:N:W, 1:N:H) = Zj;
    
% decode the encoded dc coefficients (with a 2Nx2N dct block if you put 2*N
% into lbt_dec)
Zi(1:N:W, 1:N:H) = lbt_dec(Zi(1:N:W, 1:N:H), N, s)*N;

Z = lbt_dec(Zi, N, s);

% clip Z values to be in 0 to 255 range
Z = Z + 128;
Z(Z<0) = 0;
Z(Z>255) = 0;

%fprintf(1, 'Inverse %i x %i DCT\n', N, N);
%C8=dct_ii(N);
%Z=colxfm(colxfm(Zi',C8')',C8');

return