/* Gaussian random number sequence generator by Walsh transform */
/* pseudorandom numbers ---> Guassian pseudorandom numbers      */
/* (c) 2025 cepstrum.co.jp                                      */

clear;

//SEED='a5a5a5a5';
//SEED='5a5a5a5a';
//SEED='0f0f0f0f';
SEED='f0f0f0f0';

grand('setgen', 'mt')                    // use Mersenne-Twister 
grand('setsd',  hex2dec(SEED));          // set seed

DATLEN =2^20;                            // data length = Walsh transform length (2^20=1048576)
HIST_DIV=100;                            // histragem parameter

g=zeros(1:DATLEN);                       // workspace of Walsh transform
g2=zeros(1:DATLEN);                      // workspace of Walsh transform
x=zeros(1:DATLEN);
y=zeros(1:DATLEN);

// generate random signal without DC offest
x=grand(1, DATLEN, 'lgi');               // uniform distribution
x=x-sum(x)/DATLEN;                       // remove DC offset
x=x/max(abs(x));                         // scaling

g=x;                                     // copy to Walsh routine workspace  

// Walsh transform
// converted from FORTRAN code in the old book published more than 40 years ago
// *** very slow under Scilab interpreter ***
mprintf("           |");
for i=1:log2(DATLEN)-1
  mprintf("_");
end
mprintf("|\nworking... ");
for i1=0:log2(DATLEN)-1
  mprintf("*");
  i2=2^(i1+1);
  for j=0:i2-1
    i3=DATLEN/(2^(i1+1));
    for k=0:i3-1
      n1=k*(2^(i1+1))+j;
      n2=2*k*(2^i1)+floor(j/2);
      n3=(2*k+1)*(2^i1)+floor(j/2);
      pn=(-1)^floor((j+1)/2);
      g2(n1+1)=g(n2+1)+pn*g(n3+1);
    end
  end
  g=g2;
end
mprintf("*");

y=g;          // copy output sequnece

// save as wave file
//ymax=1.05*max(abs(y));
ymax=3300;
savewave('out.wav', 0.9*y/ymax, 8000);

// plot waveform 
scf(0);
clf;
plot2d(y(1:1000), axesflag=1, style=2, rect=[0, -ymax, 1000, ymax]);

// calculate histgram 
hist_range=2.0*ymax/(2*HIST_DIV+1);
hist=zeros(1:2*HIST_DIV+1);
for i=1:DATLEN
  ix=floor((y(i)+ymax)/hist_range)+1;
  hist(ix)=hist(ix)+1;
end

// plot histgram 
scf(1);
//clf;           // clear screen
gcf().axes_size = [800, 600];
plot2d(((1:2*HIST_DIV+1)-HIST_DIV-1)*hist_range, hist, axesflag=1, rect=[-ymax, 0, ymax, 28000], style=2);
plot2d([0, 0], [0, 28000],style=35);
title fontsize 2
title("cepstrum.co.jp");
xlabel "output histgram (Gaussian) "




