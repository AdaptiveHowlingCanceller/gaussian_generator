/* Gaussian random number sequence generator by FFT */
/* (c) 2025 cepstrum.co.jp                          */

clear;                                             // all clear
DATLEN=2^20;                                       // length=2^20=1048576
rand('uniform');                                   // set uniform random number mode 
rand('seed', hex2dec('a5a50f0f'));                 // set random number seed

//==========================================================================================
x=(rand(1:DATLEN)-0.5)+%i*(rand(1:DATLEN)-0.5);    // generate complex uniform random number
y=fft(x, -1);                                      // FFT
gauss_rand=[real(y), imag(y)];                     // complex ---> real conversion
//==========================================================================================

scf(0); clf;
plot2d(gauss_rand(1:3000), axesflag=1);

scf(1); clf;
histplot(-1500:50:1500, gauss_rand, axesflag=1);

savewave('out.wav', 0.9*gauss_rand/max(abs(gauss_rand)));


