function [vlc bits huffval q_ratio] = jpegenclbtenergy(X, qstep, N, M, rise1, s, opthuff, dcbits)
    
% Same as before but using the image to work out q ratios which are
% inversely proportional to the energies of the sub images 
% JPEGENC Encode an image to a (simplified) JPEG bit stream
%
%  [vlc bits huffval] = jpegenc(X, qstep, N, M, opthuff, dcbits) Encodes the
%  image in X to generate the variable length bit stream in vlc.
%
%  X is the input greyscale image
%  qstep is the quantisation step to use in encoding
%  N is the width of the DCT block (defaults to 8)
%  M is the width of each block to be coded (defaults to N). Must be an
%  integer multiple of N - if it is larger, individual blocks are
%  regrouped.
%  rise1 is half zero step size/normal step size. Default value is 1, which
%  means the zero step is double the normal step size.
%  s is the step size for the LBT. This defaults to sqrt(2).
%  if opthuff is true (defaults to false), the Huffman table is optimised
%  based on the data in X
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
error(nargchk(2, 8, nargin, 'struct')); % the 8 was originally 6, don't know if this actually needed changing
%if ((nargout~=1) && (nargout~=3)) error('Must have one or three output arguments'); end
if (nargin<8)
  dcbits = 16; % updated the default to 16 from 8
  if (nargin<7)
    opthuff = false;
    if (nargin<6)
        rise1 = 1; %default step size 1
        if (nargin<5)
            s = sqrt(2);
            if(nargin<4)
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
 %if ((opthuff==true) && (nargout==1)) error('Must output bits and huffval if optimising huffman tables'); end
 
% DCT on input image X.
%fprintf(1, 'Forward %i x %i DCT\n', N, N);
%C8=dct_ii(N);
%Y=colxfm(colxfm(X,C8)',C8)'; 

% LBT on the input image X.
%size of the image
I = size(X, 1);
% pre filter
[Pf, ~] = pot_ii(N,s); 

t = (1+N/2):(I-N/2);% N is the DCT size, I is the image size
Xp = X; 
Xp(t,:) = colxfm(Xp(t,:),Pf);
Xp(:,t) = colxfm(Xp(:,t)',Pf)';

CN = dct_ii(N);
%Giving Xp zero mean
Xpzero = Xp-128;
Y=colxfm(colxfm(Xpzero,CN)',CN)';

% finding the q_ratios
en = dct_energies(Y,N); % energy of each sub image
en = en/sum(en(:)); % normalized by the total energy
q_ratio = 1./sqrt(en); % q ratios are inversely proportional to the square root of the energies

q_ratio = q_ratio./q_ratio(N,N); % normalising the q ratios
% Quantise to integers.
fprintf(1, 'Quantising to step size of %i\n', qstep); 
Yq=quantratio1(Y, qstep, q_ratio, N, rise1);
%Yq=quant1(Y,qstep,qstep*rise1);

% Generate zig-zag scan of AC coefs.
scan = diagscan(M);

% On the first pass use default huffman tables.
disp('Generating huffcode and ehuf using default tables')
[dbits, dhuffval] = huffdflt(1);  % Default tables.
[huffcode, ehuf] = huffgen(dbits, dhuffval);

% Generate run/ampl values and code them into vlc(:,1:2).
% Also generate a histogram of code symbols.
disp('Coding rows')
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
  fprintf(1,'Bits for coded image = %d\n', sum(vlc(:,2)));
  return;
end

% Design custom huffman tables.
disp('Generating huffcode and ehuf using custom tables')
[dbits, dhuffval] = huffdes(huffhist);
[huffcode, ehuf] = huffgen(dbits, dhuffval);

% Generate run/ampl values and code them into vlc(:,1:2).
% Also generate a histogram of code symbols.
disp('Coding rows (second pass)')
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
fprintf(1,'Bits for coded image = %d\n', sum(vlc(:,2)))
fprintf(1,'Bits for huffman table = %d\n', (16+max(size(dhuffval)))*8)

if (nargout>1)
  bits = dbits;
  huffval = dhuffval';
end

return
