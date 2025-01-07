function ret = smooth_3rd_oct( xft )
%Smoothing very similar to what Abhi obtained - only minor departures
%at higher frequencies
%x is the fft of the signal (impulse response)

%xft2=xft.^2;
xft2=xft;
lx = length(xft);
smoothed=zeros([lx 1]);
ind=0:lx-1;
i_min = round(ind/2^(1/6));
i_max = min(round(ind*2^(1/6)),lx-1);
for i=ind
   smoothed(i+1) = mean(xft2((i_min(i+1):i_max(i+1))+1));
end
ret=smoothed;
   
