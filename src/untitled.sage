var("k,y");
f(u,v,T) = sum(exp(-k*(k-1)*T/2) *(2*k-1)*(-1)^(k-v) / (factorial(v) * factorial(k-v)*(v+k-1))*exp(sum(ln((v+y)*(u-y)/(u+y)),y,0,k-1)),k,v,u)

T=0.2;
pgt1 = 1/3 * f(2,2,T).expand_sum() * 1/3 * f(3,2,T).expand_sum() * f(2,2,T).expand_sum() * f(2,2,T).expand_sum() * 2/6*1/3 * f(4,1,T).expand_sum();
pgt2 = 1/6 * 1/3 * f(2,2,T).expand_sum() * 1/3 * f(3,2,T).expand_sum() * f(2,2,T).expand_sum() * f(2,2,T).expand_sum() * 2/6*1/3*f(4,2,T).expand_sum()

print(pgt1)
print(pgt2)

T5 = T2 = T3 = T4 = 1
T1 = 0.01

pgt1 = 1/3 * f(2,2,T1).expand_sum().n() * 1/3 * f(3,2,T2).expand_sum().n() * f(2,2,T3).expand_sum().n() * f(2,2,T4).expand_sum().n() * 2/6*1/3 * f(4,1,T5).expand_sum().n();
pgt2 = 1/6 * 1/3 * f(2,2,T1).expand_sum().n() * 1/3 * f(3,2,T2).expand_sum().n() * f(2,2,T3).expand_sum().n() * f(2,2,T4).expand_sum().n() * 2/6*1/3*f(4,2,T5).expand_sum().n()

print(pgt1)
print(pgt2)
