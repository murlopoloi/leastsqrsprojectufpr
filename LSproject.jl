#Done in Julia 1.1.1, 05/20/2020 (MM/DD/YYYY)
using Plots
vx = [3 6 9 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60 63 66 69 72 75 78 81 84 87 90 93 96 99 102 105 108 111 114 117 120 123 126 129 132 135 138 141 144 147 150 153 156 159 162 165 168 171 174 177 180 183 186 189 192]
#= Data used for the project =#
#Veganism vy = [13 10 9 6 9 7 7 8 5 8 8 8 4 6 6 7 8 5 8 8 7 5 6 7 9 6 8 9 8 8 8 11 8 10 10 11 10 11 15 17 19 19 22 26 24 30 37 37 42 41 48 51 57 58 65 71 73 67 69 79 100 86 92 98]
#Vegetarianism vy = [41 35 34 27 31 19 33 25 26 35 24 29 25 24 27 24 28 24 33 30 25 22 22 22 20 22 21 21 22 18 18 18 18 18 18 17 20 18 20 22 25 19 23 27 22 24 29 30 29 28 27 33 35 28 33 35 29 28 28 34 37 29 34 38]

function lsqm(vx,vy) #Linear least square
    m = size(vx,2)
    sumyi, sumxi, sumxiyi, sumxisq, sumyisq, sumerrabs = 0, 0, 0, 0, 0, 0
    for i = 1:m
        sumyi = sumyi + vy[i]
        sumxi = sumxi + vx[i]
        sumxiyi = sumxiyi + vx[i] * vy[i]
        sumxisq = sumxisq + vx[i] ^ 2
        sumyisq = sumyisq + vy[i] ^ 2
    end
    α = m * sumxisq - sumxi ^ 2  
    if α != 0
        a₀ = 1/α * (sumyi * sumxisq - sumxi * sumxiyi)
        a₁ = 1/α * (m * sumxiyi - sumxi * sumyi)
    end
    for i = 1:m
        sumerrabs = sumerrabs + ((a₀ + a₁ * vx[i]) - vy[i]) ^ 2
    end 
    return a₀, a₁, sqrt(sumerrabs), sqrt(sumerrabs/sumyisq)
end

function expsqm(vx,vy) #Exponential least square
    m = size(vx,2)
    sumyi, sumxi, sumxiyi, sumxisq, sumyisq, sumerrabs = 0, 0, 0, 0, 0, 0
    for i = 1:m
        sumyi = sumyi + log(vy[i]) 
        sumxi = sumxi + vx[i]
        sumxiyi = sumxiyi + vx[i] * log(vy[i])
        sumxisq = sumxisq + vx[i] ^ 2
    end
    α = m * sumxisq - sumxi ^ 2  
    if α != 0
        a₀ = 1/α * (sumyi * sumxisq - sumxi * sumxiyi)
        a₁ = 1/α * (m * sumxiyi - sumxi * sumyi)
    end
    for i = 1:m
        sumerrabs = sumerrabs + (ℯ ^ a₀ * ℯ ^ (a₁ * vx[i]) - vy[i]) ^ 2 
    end 
    return ℯ^a₀, a₁, sqrt(sumerrabs), sqrt(sumerrabs/sumyisq)
end

function quadsqm(vx,vy) #Quadratic least square
    m = size(vx,2)
    sumyi, sumxi, sumxiyi, sumxisq, sumyisq, sumxicb, sumxifth, sumxisqyi, sumerrabs = 0, 0, 0, 0, 0, 0, 0, 0, 0
    for i = 1:m
        sumyi = sumyi + vy[i] 
        sumxi = sumxi + vx[i]
        sumxiyi = sumxiyi + vx[i] * vy[i]
        sumxisq = sumxisq + vx[i] ^ 2
        sumxicb = sumxicb + vx[i] ^ 3
        sumxifth = sumxifth + vx[i] ^ 4
        sumxisqyi = sumxisqyi + vx[i] ^ 2 * vy[i]
    end
    A = (m * sumxiyi - (sumxi * sumyi))/(m * sumxisq - (sumxi) ^ 2)
    B = (m * sumxicb - (sumxisq * sumxi))/(m * sumxisq - (sumxi) ^ 2)
    C = sumxicb - 1/m * (sumxi * sumxisq)
    a₂ = (m * sumxisqyi - (sumxisq * sumyi) - m*A*C)/(m * sumxifth - (sumxisq) ^ 2 - m * B * C)  
    a₀ = 1/m * (sumyi - A * sumxi - a₂ * (sumxisq - B * sumxi))
    a₁ = A - B * a₂
    for i = 1:m
        sumerrabs = sumerrabs + (a₀ + a₁ * vx[i] + a₂ * vx[i] ^ 2 - vy[i]) ^ 2 
    end
    return a₀, a₁, a₂, sqrt(sumerrabs), sqrt(sumerrabs/sumyisq)
end
    
function sqexpsqm(vx,vy) #Quadratic exponential least square
    m = size(vx,2)
    sumyi, sumxi, sumxiyi, sumxisq, sumyisq, sumxicb, sumxifth, sumxisqyi, sumerrabs = 0, 0, 0, 0, 0, 0, 0, 0, 0
    for i = 1:m
        sumyi = sumyi + log(vy[i]) 
        sumxi = sumxi + vx[i]
        sumxiyi = sumxiyi + vx[i] * log(vy[i])
        sumxisq = sumxisq + vx[i] ^ 2
        sumxicb = sumxicb + vx[i] ^ 3 
        sumxifth = sumxifth + vx[i] ^ 4
        sumxisqyi = sumxisqyi + vx[i] ^ 2 * log(vy[i])
    end
    A = (m * sumxiyi - (sumxi * sumyi))/(m * sumxisq - (sumxi) ^ 2)
    B = (m * sumxicb - (sumxisq * sumxi))/(m * sumxisq - (sumxi) ^ 2)
    C = sumxicb - 1/m * (sumxi * sumxisq)
    a₂ = (m * sumxisqyi - (sumxisq * sumyi) - m*A*C)/(m * sumxifth - (sumxisq) ^ 2 - m * B * C)  
    a₀ = 1/m * (sumyi - A * sumxi - a₂ * (sumxisq - B * sumxi))
    a₁ = A - B * a₂
    for i = 1:m
        sumerrabs = sumerrabs + (ℯ ^ (a₀ + a₁ * vx[i] + a₂ * vx[i] ^ 2) - vy[i]) ^ 2 
    end
    return a₀, a₁, a₂, sqrt(sumerrabs), sqrt(sumerrabs/sumyisq)
end
