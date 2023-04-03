%Develloped by:Khashayar Namdar, reviewed by Seyedashkan Miraftabi
clear;
clc;
close all;
N=64; %Number of subcarries, it should be power of 2. 64 and more is preffered.
M=4; %Number of symbol states, M=4 for QPS (it should be power of 2 and less than N)
L=1; %L factor
paprdb1s=1;
paprdb1=100;
paprdb=0;


while paprdb1>paprdb-4
                                   % Normal OFDM
%QPSK modulation
r=floor(M*rand(N,1)); %floor: Round toward negative infinity
%rand:Normally destributed between 0and1, rand(N,1)=N by 1 matrix
%r[N,1] and elements are 0,1,2,3 

bexp = pskmod(r,M);
%y = pskmod(x,M) modulates the input signal, x, using phase shift keying (PSK)
%with modulation order M.
%bex[N,1] and elements are (-1,-1),(-1,1),(1,-1),(1,1) in complex numbers

% Calculating |IFFT| square and displaying |IFFT|
ibexp=ifft(bexp);%X = ifft(Y) computes the inverse discrete Fourier 
%transform of Y using a fast Fourier transform algorithm. X is the same size as Y.

mibexp=abs(ibexp);%(a^2+b^2)^0.5 (Complex Magnitude)
%mibexp[N,1], real numberbetween 0,1

smibexp=mibexp.^2; %smibexp[N,1] square complex magnitude

%Calculation of PAPR
papr=(max(smibexp))/(mean(smibexp)); %papr formulation

paprdb=10*log(papr); %db formulation
                                  % OFDM modified by Selective Mapping
                                  % Technique

% Normalized Riemann Matrix used for SLM technique(Selected Mapping)                                
rm=gallery('riemann',N); %rm[N,N]
b=rm/N; %b[N,N]
%riemann   Matrix associated with the Riemann hypothesis
%A = gallery('riemann',n) returns an n-by-n matrix for which the Riemann hypothesis is true if and only if
%det(A)=O(n!n??1/2+?)
%for every ? > 0.
%The Riemann matrix is defined by:
%A = B(2:n+1,2:n+1) 

% Elementwise Multiplication(Each element of each rown of b is multipled by bexp corresponding element)
for i=1:N
    for j=1:N
    bexp1(i,j)=b(i,j).*bexp(j,1); %bexp1[N,N]
    end;
end;

% Calculating |IFFT| square
ibexp1=ifft(bexp1);
mibexp1=abs(ibexp1);
smibexp1=mibexp1.^2;

% Calculating PAPR for each of the N blocks
for i=1:N
papr1(i,1)=(max(smibexp1(i,:)))/(mean(smibexp1(i,:)));
%papr1(N,1)
end;

% Finding the block with minimum PAPR
xm=1;
for i=2:N
    if (papr1(i,1)<papr1(xm,1))
     xm=i;
    end;
end;

% Minimum PAPR in dB
paprdb1=10*log(papr1(xm,1));

a=bexp;
LN=floor(L*N);
at=a';
aa=[at(1:N) zeros(1,LN-N)]'; %Zero Padding if L=1, then aa[N,1]


x=ifft(aa);  %Generating OFDM signal and calculating PAPR
x_mag=abs(x); %x_mag[NL,1]
 
paprs=max(x_mag.^2)/mean(x_mag.^2);
paprdbs=10*log(paprs);
x_max=0.7*max(x_mag);%Clipping Level=0.7Max(Clipping Just the Maximum number)

for j=1:LN                                   %Clipping the signals above threshold(here 0.2)
    if(x_mag(j,1)>x_max)
        x_mag1(j,1)=x_max;
    else
        x_mag1(j,1)=x_mag(j,1);
    end;    
end;
%x_mag1[LN,1] and it is clipped

subplot(3,1,2),stem(x_mag1);
%subplot(m,n,p) divides the current figure into an m-by-n grid 
%and creates axes in the position specified by p
%stem(Y) plots the data sequence, Y, as stems that extend from a baseline along the x-axis. 
%The data values are indicated by circles terminating each stem.
%If Y is a vector, then the x-axis scale ranges from 1 to length(Y)
xlim([0 LN]); 
%xlim(limits) sets the x-axis limits for the current axes or chart


%Filtering the clipped signal 
h=[ones(1,N) zeros(1,LN-N)]';%h=[1 1 1 1 1 1 0 0 0 ](N ones +zero padding)
x_mag2=conv(x_mag1,h);
%w = conv(u,v) returns the convolution of vectors u and v

subplot(3,1,3),stem(x_mag2);
xlim([0 2*LN]); 

                                           %Calculating PAPR of Clipped and Filtered signal
papr1s=max(x_mag2.^2)/mean(x_mag2.^2);
paprdb1s=10*log(papr1s);

end;

disp('PAPR of normal OFDM=');
disp(paprdb);
disp('PAPR of SLM modified OFDM=');
disp(paprdb1);
disp('PAPR of clipped OFDM=');
disp(paprdb1s);

subplot(2,1,1),stem(mibexp),title('Normal OFDM signal');
xlim([0 N]); 
subplot(2,1,2),stem(mibexp1(xm,:)),title('SLM modified OFDM signal');
xlim([0 N]); 
figure,subplot(2,1,1),stem(x_mag),title('Normal OFDM signal with zero padding');
xlim([0 LN]); 
subplot(2,1,2),stem(x_mag1),title('Clipped OFDM signal');
xlim([0 LN]); 
eff1=(1-paprdb1/paprdb)*100;
eff2=(1-paprdb1s/paprdb)*100;
disp('%SLM Efficiency');
disp(eff1);
disp('%Clipping+Filtering Efficiency');
disp(eff2);
if(eff1>eff2)
    disp('Then SLM is more effective for PAPR reduction');
else
    
    disp('Then Clipping+Filtering is more effective for PAPR reduction');
end;