x = 5.5

sum = 0.0;
for z = 0:20
    one = x^z * gamma(z+1) / (factorial(z) * (1+x)^(z+1))
    two = poisspdf(z,x)
end

sum