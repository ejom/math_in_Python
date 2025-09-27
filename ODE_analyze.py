from scipy.integrate import solve_ivp
from scipy.special import gamma, airy
import numpy as np
import matplotlib.pyplot as plt

eq1 = "y'' + y' + y"
eq2 = "y''-ty"
eq3 = "(t^2+1)y'''+ty''+t(y')^2+y"

#we need a function that turns the string equations into a function of the form y' = f(t, y) for our solve_ivp to use

def f1(t, y):
    return [y[1], -y[1] - y[0]]
y1_ini = [1, 0]
grad1 = [[0, 1], [-1, -1]]

#y(x) = 1/3 e^(-x/2) (sqrt(3) sin((sqrt(3) x)/2) + 3 cos((sqrt(3) x)/2))
def explSol1(t):
    return (1/3)*np.exp(-t/2)*(np.sqrt(3)*np.sin((np.sqrt(3)*t)/2)+3*np.cos((np.sqrt(3)*t)/2))

def f2(t, y):
    return [y[1], t*y[0]]
y2_ini = [1 / 3**(2/3) / gamma(2/3), -1 / 3**(1/3) / gamma(1/3)]
grad2 = lambda t, y: [[0, 1], [t, 0]]

def f3(t, y):
    return [y[1], y[2], (-t*y[2]-t*(y[1])**2-y[0])/(t**2+1)]
y3_ini = [1, 0, 0]
grad3 = lambda t, y: [[0, 1, 0], [0, 0, 1], [-1/(t**2+1), -t*y[1]/(t**2+1), -t/(t**2+1)]]

t_span = [0, 10]

sol1 = solve_ivp(f1, t_span, y1_ini, method='Radau', jac=grad1, rtol=1e-8, atol=1e-8)
sol2 = solve_ivp(f2, t_span, y2_ini, method='Radau', jac=grad2, rtol=1e-8, atol=1e-8)
sol3 = solve_ivp(f3, t_span, y3_ini, method='Radau', jac=grad3, rtol=1e-8, atol=1e-8)

# Create separate figures for each ODE solution

# Figure 1: Damped oscillator (y'' + y' + y = 0)
plt.figure()
plt.subplot(1, 2, 1)
plt.plot(sol1.t, sol1.y[0])
plt.xlabel('t')
plt.ylabel('y')
plt.title('Numerical solution: y\'\' + y\' + y = 0')

plt.subplot(1, 2, 2)
plt.plot(sol1.t, explSol1(sol1.t))
plt.xlabel('t')
plt.ylabel('y')
plt.title('Explicit solution: y\'\' + y\' + y = 0')

print("Sol1 max error")
sol1_err = np.max(np.abs(sol1.y[0] - explSol1(sol1.t)))
print(sol1_err)
print(np.max(np.abs(sol1.y[0])))

# Figure 2: Airy equation (y'' - t*y = 0)
plt.figure()
plt.subplot(1, 2, 1)
plt.plot(sol2.t, sol2.y[0], label='y(t) Numerical')
plt.xlabel('t')
plt.ylabel('y')
plt.title('Numerical solution: y\'\' - t*y = 0 (Airy Equation)')

# Compare with analytical Airy solution
plt.subplot(1, 2, 2)
plt.plot(sol2.t, airy(sol2.t)[0], label='Airy function solution')
plt.xlabel('t')
plt.ylabel('Airy function')
plt.title('Comparison with Airy Function')

print("Sol2 max error")
sol2_err = np.max(np.abs(sol2.y[0] - airy(sol2.t)[0]))
print(sol2_err)
print(np.max(np.abs(sol2.y[0])))

# Figure 3: Third-order equation (y''' + t*y'' - t*(y')^2 - y = 0)
plt.figure()
plt.plot(sol3.t, sol3.y[0], label='y(t)')
plt.xlabel('t')
plt.ylabel('y')
plt.title('Numerical Solution: (t^2+1)y\'\'\'+ty\'\'+t(y\')^2+y=0')
plt.ylim([-100, 2])

plt.show()



