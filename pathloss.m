function [dB] = pathloss(d,fc)
dB=(36.7*log10(d)+22.7+26*log10(fc));
end

