function nbits = jpegbits(vlc, opthuff, dwt)
%JPEGBITS Calculates number of bits of a huffman code
%   Detailed explanation goes here
if opthuff
    nbits = vlctest(vlc) + 1424;
else
    nbits = vlctest(vlc);
end

if dwt
    nbits = nbits + 88;
else
    nbits = nbits + 120;
end

end
