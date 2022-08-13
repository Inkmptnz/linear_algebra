function while_loop(n)
    i = 0
    s = 0
    while i < n
        s+=i
        i += 1
    end
    s
end

function for_loop(n)
    s = 0
    for i in 0:(n-1)
        s += i
    end
    s
end

function sumi(n)
    return (n * (n - 1)) / 2
end