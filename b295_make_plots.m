mix_resistors = 1E3 * [33 39 33 33 33 33 33 39 33 33];

%There is a slight loss from taking the inputs to the mixer from
%the 220 ohm current limiting resistors at the outputs of the filter
%op amps instead of the op amps directly, but you still get like
%99 percent or 99.1 percent of the signal, so I decided not to bother.
%I left these lines here in case you want to put it back in.
%par50_mix = 50E3*mix_resistors ./ (50E3+mix_resistors);
%postatten = par50_mix ./ (200+par50_mix);

preatten = [480/3780 330/3780 330/3780 330/3780 
            330/3780 330/3780 330/3780 480/3780];
postamp = 56E3 ./ mix_resistors;

% Don left a note on the schematic to trim the 3.5 kHz
% section at 3.2 kHz, so we will treat it as a 3.2 kHz 
% section; maybe Don changed his mind about the filterbank
% design after he already had the front panels
% fabricated?
targets = [200 350 500 700 1000 1400 2000 3200];

c1arr = [0.047E-6 0.022E-6 0.022E-6 0.01E-6 ...
         0.01E-6 0.0047E-6 0.0047E-6 0.0022E-6];
c2arr = [0.01E-6 0.0047E-6 0.0022E-6 0.0015E-6 910E-12 910E-12 470E-12 470E-12];
c3arr = c1arr;

r1arr = 15E3 * ones(1,8);
% For inital guess, calculate R2 based on R3=R1 assumption
r2guess = 1 ./ ((2*pi*targets) .^ 2 .* r1arr .* c1arr .* c2arr);

% For the next two code blocks, comment in or out the
% "initial guess" and "calibrated" lines as needed depending
% on which plot you are trying to make

% inital guess for R2
% r2arr = [89.8 133 140 230 186 201 191 159] * 1E3;
% calibrated R2 values
r2arr = [93 136.6 149.5 242 199.5 209.8 200.5 171.2] * 1E3;

% inital guess for R3
% r3arr = [15 15 15 15 15 15 15 15] * 1E3;
% calibrated R3 values
r3arr = [13.98 14.04 13.15 12.41 12.75 12.74 13.4 11.3] * 1E3;
%r3arr = r1arr; % to double check formulas R1=R3 case

f = 50:1:10000;
s = 2*pi*j*f;

h = zeros(10,length(f));;

% lowpass 3rd-order Sallen-Key
r = 56E3;
c1 = 0.0015E-6; c2 = 0.22E-6; c3 = 0.047E-6; 
h(1,:) = (39/33) ./ ...
  (r^3*c1*c2*c3*s .^ 3 + 2*r^2*c1*(c2+c3)*s .^ 2 ...
   + r*(3*c1 + c3)*s + 1);

% highpass 3rd-order Sallen-Key
r1 = 220E3; r2 = 1.5E3; r3 = 6.8E3; 
c = 0.0022E-6;
h(10,:) = (39/33) * r1*r2*r3*c^3*s .^ 3 ./ ...
  (r1*r2*r3*c^3*s .^ 3 + r2*(r1+3*r3)*c^2*s .^ 2 ...
   + 2*(r2+r3)*c*s + 1);

% middle bands
for k = 1:8,
  c1 = c1arr(k); c2 = c2arr(k); c3 = c3arr(k);
  r1 = r1arr(k); r2 = r2arr(k); r3 = r3arr(k);
  b2 = r2 * r3 * c1 * c2 + r1 * r3 * c1 * (c2 + c3);
  b1 = (r2 + r3) * c2 + r3 * c1;
  a3 = r1 * r2 * r3 * c1 * c2 * c3;
  a2 = r1 * (r2 + r3) * c2 * c3;
  a1 = r1 * (c2 + c3);

  numer = b2 * s .^ 2 + b1 .* s;
  denom = a3 * s .^ 3 + a2 * s .^2 + a1 * s + 1;
  zer = roots([b2 b1])
  zer_emp_hz(k+1) = zer/(2*pi); 
  % predicted zero
  zer_pred = -1/(r1*c1)
  pol = (roots([a3 a2 a1 1]))
  real_pol_emp_hz = real(pol(1))/(2*pi);
  % empirical natural frequency from computed poles
  fnat_emp(k+1) = sqrt(pol(2)*pol(3))/(2*pi);
  wcp = sqrt(pol(2)*pol(3));
  q_emp(k+1) = -wcp / real(pol(2) + pol(3));
  h(k+1,:) = numer ./ denom;
  a_emp(k+1) = max(abs(h(k+1,:)));
  h(k+1,:) = (39/33) * preatten(k) * h(k+1,:);
  fnat_pred(k+1) = 1 / (2*pi*(c1*c2*r1*r2).^(1/2));
  q_pred(k+1) = sqrt(r2*c1/(r1*c2));
  a_pred(k+1) = r2/r1 + (1+c1/c2);
end

for k=1:10,
  hpa(k,:) = postamp(k) *  h(k,:);
end

gainsdb = 20*log10(abs(hpa));
[gains,maxindex] = max(gainsdb,[],2);
peakfreqs = f(maxindex);
% display for calibration
disp(peakfreqs(2:9));
disp(gains(2:9)');

samesum = sum(hpa,1);
flipsum = sum(hpa([1 2 4 6 8],:),1) - sum(hpa([3 5 7 9 10],:),1);
ahpa = abs(hpa)';
semilogx(f,max(gainsdb,-10)); axis tight; grid
xlabel('Frequency in Hz')
ylabel('Magnitude in dB')

% Use 
% print -deps initial_guess_freqresp.eps
% or
% print -deps -color calibrated_freqresp.eps
% to save image as needed

% Summed output plot
%semilogx(f,20*log10([samesum;flipsum]'),["-","."]'); axis tight; grid
%xlabel('Frequency in Hz')
%ylabel('Magnitude in dB')
%print -deps -color summed_output.eps