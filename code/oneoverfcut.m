function y = oneoverfcut(a, r, c, cut)
%
% ONEOVERF - returns an array of 1/f noise
%
% usage:
%	ONEOVERF(a, n) returns a square 2D array
%	of noise with a 1/(f^a) amplitude spectum
%	ONEOVERF(a, m, n) returns an m-by-n array
%
% see also: spike, oneoverf3d, spike3d
%
% Lawrence K. Cormack

% history:
% 11/12/2002    lkc     There seems to be an orientation bias when the array size is odd.
%                       The called function spike() is centering its peak correctly vis-a-vis 
%                       the shifted Fourier transform, yet the bias remains... Perhaps this is
%                       an inherent short coming of the fft algorithm.
%                       at any rate, I'm adding a cludge to pad up to the next even pixel and 
%                       then trim.
% 11/12/2000    lkc     cleaned it up 
% 3/4/01        lkc     added variable exponent
% 4/1/02        lkc     exponent is now handled in spike()
% 4/1/02        lkc     changed the way the spectrum is generated to the more intuitive
%                       amp/phase based-method, but it makes no functional difference.
% 2/19/10       tbc     changed "i" to "1i" per matlab speed suggestion
 
%*** deal with oddity ***
codd = 0; rodd = 0;
if nargin==2 c=r; end
if isodd(c)
	c=c+1;
	codd = 1;
end
if isodd(r)
	r=r+1;
	rodd = 1;
end

% ******** make the noise ************************
skirt = spike(a, r, c);	                            % the amp. spectrum
skirt(length(skirt)/2-cut:length(skirt)/2+cut, length(skirt)/2-cut:length(skirt)/2+cut) = 0;
y = skirt.*exp(1i.*2.*pi.*rand(r, c));               % noise on the complex frequency plane
y=normalize(real(ifft2(fftshift(y))));              % and over to the space domain.

%*** deal with oddity again ***
if codd y=trimr(y, 1); end
if rodd y=trimb(y, 1); end

%*********************************************** the end