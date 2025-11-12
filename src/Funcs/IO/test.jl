
# function that builds a full expanded string for an affine expression
expr=model[:CV2magexporg][7,50,31]
expr1=model[:CV2magexporg][7,20,31]


keysss=collect(intersect(keys(expr.terms), keys(expr1.terms)))
sf=length(keysss)
a=zeros(sf)
b=zeros(sf)
for j=1:sf
    i=keysss[j]
    a[j]=expr.terms[i]
    b[j]=expr1.terms[i]
end

Plots.plot(a)
Plots.plot!(b)

i=collect(keys(expr.terms))[1]
j=1

findmax(b)
keysss[115]

raw["QActBus"][21]

intersect(keys(expr.terms), keys(expr1.terms))