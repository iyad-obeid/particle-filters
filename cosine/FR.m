g = [13, 6.9, 11, 1, 4, 16, 4, 16, 7.6, 11.5]';

N = [391, 151, 164, 25, 58, 491, 127, 774, 189, 198]';

[a,b]=sort(g);
N=N(b);

a=2*a;

plot(a,700*N/30,'x');