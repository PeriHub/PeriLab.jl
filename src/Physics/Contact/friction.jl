using BenchmarkTools

x= zeros(5000)
t=0
function mod(x,y,z,t)
    return 2
end

@btime for i in 1:5000
    E.*x
end

@btime  for i in 1:5000
    mod(x,y,z,t).*x
end