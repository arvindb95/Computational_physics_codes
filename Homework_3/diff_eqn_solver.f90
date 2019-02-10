function f(x,y) result(z)
        double precision x, y, z
        z = -x*y
end function f

function euler_method(x0,y0,x,stepsize) result(e)
        double precision f, x0, y0, x, stepsize, x1, y1, e
        integer n, i
        n = (x - x0)/stepsize
        do 10 i = 0, n
                y1 = y0 + stepsize*f(x0, y0)
                x1 = x0 + stepsize
                x0 = x1
                y0 = y1
        10 continue
        e = y1
end function euler_method

function partialx(x0,y0,stepsize) result(ddx)
        double precision f, x0, y0, stepsize, ddx
        ddx = (f(x0-2.0*stepsize,y0)-8.0*f(x0-stepsize,y0)+8.0*f(x0 + stepsize,y0)-f(x0+2.0*stepsize,y0))/(12.0*stepsize)
end function partialx

function partialy(x0,y0,stepsize) result(ddy)
        double precision f, x0, y0, stepsize, ddy
        ddy = (f(x0,y0-2.0*stepsize)-8.0*f(x0,y0-stepsize)+8.0*f(x0,y0+stepsize)-f(x0,y0+2.0*stepsize))/(12.0*stepsize)
end function partialy

function taylor_method(x0,y0,x,stepsize) result(t)
        double precision f, x0, y0, x, stepsize, x1, y1, t, partialx, partialy, h
        integer n, i
        h = 10E-7
        n = (x - x0)/stepsize
        do 20 i = 0, n
                y1 = y0 + stepsize*f(x0,y0) + (stepsize**2.0)*(partialx(x0,y0,h) + f(x0,y0)*partialy(x0,y0,h))
                x1 = x0 + stepsize
                x0 = x1
                y0 = y1
        20 continue
        t = y1
end function taylor_method

function adam_bashford_method(x0,y0,x,stepsize) result(a)
        double precision f, x0, y0, x, stepsize, x1, y1, x2, y2, a
        integer n, i
        n = (x - x0)/stepsize
        x1 = x0 + stepsize
        y1 = y0 + stepsize*f(x0,y0)
        !print*,x1
        !print*,y1
        do 30 i = 1 ,n
                y2 = y1 + stepsize*(1.5*f(x1,y1) - 0.5*f(x0,y0))
                x2 = x1 + stepsize
                x0 = x1
                x1 = x2
                y0 = y1
                y1 = y2
        30 continue
        a = y2
end function adam_bashford_method

function runge_kutta(x0,y0,x,stepsize) result(r)
        double precision f, x0, y0, x, x1, y1, stepsize, k, r
        integer n, i
        n = (x - x0)/stepsize
        do 40 i = 1,n
                k = stepsize*f(x0,y0)
                y1 = y0 + stepsize*f(x0 + stepsize*0.5, y0 + k*0.5)
                x1 = x0 + stepsize
                x0 = x1
                y0 = y1
        40 continue
        r = y1
end function runge_kutta

program solve_diff_eqn

implicit none
double precision x0,y0,x,stepsize, y_real_1, e, t, a, r, y_real_3
double precision euler_method, taylor_method, adam_bashford_method, runge_kutta 

y_real_1 = EXP(-0.5)
y_real_3 = EXP(-4.5)

x0 = 0
y0 = 1
x = 3
stepsize = 0.001

e = euler_method(x0, y0, x, stepsize)

print*, e
print*, e - y_real_3
print*,"###########################################"

x0 = 0
y0 = 1
t = taylor_method(x0, y0, x, stepsize)

print*,t
print*,t - y_real_3
print*,"###########################################"

x0 = 0
y0 = 1
a = adam_bashford_method(x0, y0, x, stepsize)

print*,a
print*,a - y_real_3
print*,"###########################################"

x0 = 0
y0 = 1
r = runge_kutta(x0, y0, x, stepsize)

print*,r
print*,r - y_real_3

end program solve_diff_eqn
