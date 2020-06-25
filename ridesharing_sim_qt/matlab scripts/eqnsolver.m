cap = 3
x = 2.5


coeffs = zeros(1,cap+2);
coeffs(1) = 1;
coeffs(2) = -(1+x)/x;
coeffs(cap+2) = 1/x;

slns = roots(coeffs)

slnsabs = abs(slns)

assert(slnsabs(1) > 1)