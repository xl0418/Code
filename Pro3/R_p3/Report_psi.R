library(POJCE)

# simulate data sets for m loc

psi = 0
sA = 1000
sD = 1000
ticks = 200000
L_size = 51
v = 0.001
seed = 100

phi = 0
result = pojce(phi=phi, psi=psi, sA=sA, sD=sD, ticks=ticks, randomSeed=seed, L=L_size, v=v)
print(tail(result$events))

phi = 0.5
result = pojce(phi=phi, psi=psi, sA=sA, sD=sD, ticks=ticks, randomSeed=seed, L=L_size, v=v)
print(tail(result$events))

phi = 1
result = pojce(phi=phi, psi=psi, sA=sA, sD=sD, ticks=ticks, randomSeed=seed, L=L_size, v=v)
print(tail(result$events))


phi = 1
psi = 0
result1 = pojce(phi=phi, psi=psi, sA=sA, sD=sD, ticks=ticks, randomSeed=seed, L=L_size, v=v)
# print(tail(result$events))

psi = 0.5
result2 = pojce(phi=phi, psi=psi, sA=sA, sD=sD, ticks=ticks, randomSeed=seed, L=L_size, v=v)
# print(tail(result$events))

psi = 1
result3 = pojce(phi=phi, psi=psi, sA=sA, sD=sD, ticks=ticks, randomSeed=seed, L=L_size, v=v)
# print(tail(result$events))


