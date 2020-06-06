
pfulls = zeros(1,65);
curi = 1;

for cap = 1:8
    for x = 0.5:0.5:cap-0.5
        cap
        x
        
        zk = zeros(1,cap);

        for k = 1:cap
            zk(k) = - cap / x * lambertw(- x / cap * exp(-x / cap) * exp(- 1i * 2*pi * k / cap));
        end

        zk;


        syms pgf
        syms z
        syms temparr

        temparr = sym(zeros(1, cap));

        for l = 1:(cap-1)
            temparr(l) = (z - zk(l)) / (1-zk(l));
        end

        temparr;

        pgf = sym(1);

        for j = 1:(cap-1)
            pgf = pgf * temparr(j);
        end

        pgf = pgf * (cap - x) * (z - 1) / (z^cap * exp(x * (1-z)) - 1);

        syms lastderiv

        lastderiv = pgf;
        numcalc = 20;
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
        sump;

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
        curi = curi+1;
    end
end

pfulls

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