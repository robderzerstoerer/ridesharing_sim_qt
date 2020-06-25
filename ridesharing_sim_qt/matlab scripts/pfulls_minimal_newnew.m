poissonbuses = true;
pfulls = zeros(1,65);
pfullsnew = zeros(1,65);
curi = 1;

for cap = 1:8
    for x = 0.5:0.5:cap-0.5
        cap
        x
        
        zk = zeros(1,cap);

        if ~poissonbuses
            for k = 1:cap
                zk(k) = - cap / x * lambertw(- x / cap * exp(-x / cap) * exp(- 1i * 2*pi * k / cap));
            end
        end
        

        syms pgf
        syms z
        syms temparr

        temparr = sym(zeros(1, cap));

        if ~poissonbuses
            for l = 1:(cap-1)
                temparr(l) = (z - zk(l)) / (1-zk(l));
            end
            
            pgf = sym(1);
            for j = 1:(cap-1)
                pgf = pgf * temparr(j);
            end

            pgf = pgf * (cap - x) * (z - 1) / (z^cap * exp(x * (1-z)) - 1);
        
        else
            coeffs = zeros(1,cap+2);
            coeffs(1) = 1;
            coeffs(2) = -(1+x)/x;
            coeffs(cap+2) = 1/x;

            slns = roots(coeffs);

            slnsabs = abs(slns);

            assert(slnsabs(1) > 1);
            zs = slns(1);
            
            pgf = (zs - 1) / (zs - z);
        end

        syms lastderiv

        lastderiv = pgf;
        numcalc = 2*cap - 1;
        sump = 0.0;
        pj = zeros(1, numcalc+1);
        
        pj(1) = vpa(subs(lastderiv, z, 0));
        sump = pj(1);

        for j = 1:numcalc
            j
            lastderiv = diff(lastderiv);
            pj(j + 1) = vpa(subs(lastderiv, z, 0)) / factorial(j);
            sump = sump + pj(j+1);
        end

        pj
        sump

        ppz = zeros(numcalc+1-cap,numcalc+1);

        for zz=0:(numcalc-cap)
           summ = 0.0;
            for m = 0:cap
                summ = summ + pj(m+1);
            end
            ppz(zz+1,zz+1) = summ;
            nmax = 0;
            if (numcalc+zz-cap) < numcalc
                nmax = numcalc+zz-cap;
            else
                nmax = numcalc;
            end
            for n=(zz+1):nmax
                ppz(zz+1,n+1) = pj(cap+(n-zz)+1);
            end
        end
        ppz;

        sumzz = 0.0;

        for zz = 1:(numcalc-cap)
            sumd = 0.0;
            for d = 1:(zz-1)
                sumd = sumd + d * ppz(zz+1,cap+d+1);
            end
            sumn = 0.0;
            for n = 0:(cap+zz-1)
                sumn = sumn + ppz(zz+1,n+1);
            end
            sumzz = sumzz + poisspdf(zz,x) * (sumd + zz*(1-sumn));
        end

        pfulls(curi) = sumzz / x;
        
        
        
        % NEWNEW CALCULATION
        numcalc_dzq = 50;
        
        
        sumtotal = 0.0;
        
        %sum 1
        for z = (cap+1):numcalc_dzq
            sumqp = 0.0;
            for qp = 0:cap
                sumqp = sumqp + pj(qp+1);
            end
            kz = 0;
            if ~poissonbuses
                kz = poisspdf(z,x);
            else
                kz = x^z * gamma(z+1) / (factorial(z) * (1+x)^(z+1));
            end
            sumtotal = sumtotal + kz * (z-cap) * sumqp;
        end
        
        %sum 2
        for qp = (cap+1):(2*cap-1)
            for z = (2*cap-qp):numcalc_dzq
                kz = 0;
                if ~poissonbuses
                    kz = poisspdf(z,x);
                else
                    kz = x^z * gamma(z+1) / (factorial(z) * (1+x)^(z+1));
                end
                sumtotal = sumtotal + kz * (qp+z-2*cap) * pj(qp+1); 
            end
        end
        
        %sum 3
        for z=1:numcalc_dzq
            sumqp = 0.0;
            for qp = 0:(2*cap-1)
                sumqp = sumqp + pj(qp+1);
            end
            kz = 0;
            if ~poissonbuses
                kz = poisspdf(z,x);
            else
                kz = x^z * gamma(z+1) / (factorial(z) * (1+x)^(z+1));
            end
            sumtotal = sumtotal + z*kz*(1-sumqp);
        end
        
        pfullsnew(curi) = sumtotal / x;
        
        curi = curi+1;
    end
end

pfulls
pfullsnew

%{
E = 0.0;
for n = 0:numcalc
    E = E + n * pj(n+1);
end
E

sumtwo = 0.0;
for n = 0:cap
    sumtwo = sumtwo + (n - cap) * pj(n+1);
end

pfullappr = 1 - (cap + sumtwo) / E;
pfullappr

sumthree = 0.0;
for n = (cap+1):numcalc
    sumthree = sumthree + (n-cap) * pj(n+1);
end

pfullappr2 = (sumthree) / E
%} 