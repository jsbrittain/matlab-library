function res = ht_freqres( epochs, rate, opt_str )

% Get epoch length
len = size(epochs,1);

% Construct sinusoid and replicate to 20 epochs (single channel)
dat = repmat( sin(2*pi*(rate/4)*[0:len-1]/rate)', 1, 20 );

% Analyse sinusoid using given parameters
[sp11,sp22,sp12,params] = ht_sp2_epochs( dat, dat, rate, opt_str);

% Evaluate response at epoch mid-point
response = 10*log10(sp11(:,round(end/2)));
response = response - max(response);

% Find -30 dB points
f1 = find(diff(response>-30)>0);
f2 = find(diff(response>-30)<0);

% Determine effective frequency resolution
df = diff(params.freqs([1 2]));
res = df*(f2 - f1);
