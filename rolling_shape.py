#!/usr/bin/python3

import numpy as np
import argparse
from scipy.interpolate import CubicSpline
from matplotlib import pyplot as plt

##########################################################################################

def ellipse(w1=1.0, w2=2.0, nthetas=10000):
    '''
    Generate the polar representation r = r(theta) of an ellipse
    w1 = length of semi axis 1
    w2 = length of semi axis 2
    nthetas = number of thetas
    '''
    thetas = np.linspace(0.0, 2.0*np.pi, nthetas, endpoint=True)
    rs = np.sqrt(w1**2*np.cos(thetas)**2 + w2*np.sin(thetas)**2)
    return thetas, rs

##########################################################################################

def polygon(n, w=1.0, nthetas=10000):
    '''
    Generate the polar representation r = r(theta) of a regular polygon
    n = number of polygon sides
    w = distance from centre to edge midpoint
    nthetas = number of thetas
    '''
    thetas = np.linspace(0.0, 2.0*np.pi, nthetas, endpoint=True)
    rs = np.zeros(nthetas)
    for j in range(nthetas):
        theta = thetas[j]
        #shift theta to the range [-pi/n, pi/n]
        t = theta - 2.0*np.pi/n*int((theta + np.pi/n)/(2.0*np.pi/n))
        rs[j] = w/np.cos(t)
    return thetas, rs

##########################################################################################

def polar_to_cartesian(theta, r, centre):
    '''
    Convert from polar coordinates (theta, r) to Cartesian coordinates (x, y)
    theta = angle measured from negative y axis
    r = distance from the origin
    centre = Cartesian coordinates of the origin
    '''
    #shift theta because angle measured from negative y axis
    t = theta - np.pi/2.0
    q = np.array([r*np.cos(t), r*np.sin(t)])
    q = q.T
    q += centre
    return q

##########################################################################################

def shift_input(f, a, b):
    '''
    Transform function f(x) so its argument is shifted to the range [a, b]
    f = function to be transformed
    a = lower limit of range
    b = upper limit of range
    '''
    def g(t):
        x = t - (b - a)*np.floor((t - a)/(b - a))
        return f(x) 
    return g

##########################################################################################

def integrate(xs, f):
    '''
    Calculate antiderivative of function f using Simpson's rule
    xs = quadrature points
    f = function to be integrated
    '''
    n = len(xs)
    ys = np.zeros(n)
    for j in range(1, n):
        dx = xs[j] - xs[j - 1]
        fa = f(xs[j - 1])
        fb = f(xs[j])
        fab2 = f(0.5*(xs[j - 1] + xs[j]))
        ys[j] = ys[j - 1] + (dx/6.0)*(fa + fb + 4.0*fab2)
    #ys = values of antiderivative at quadrature points
    return ys

##########################################################################################

def gudermannian(t):
    '''
    Gudermannian function
    t = argument
    '''
    g = np.arctan(np.sinh(t))
    return g

##########################################################################################

def test_quadrature():
    '''
    Test whether the quadrature agrees with the exact solution for a square
    '''
    w = 2.0
    thetas, rs = polygon(4, w)
    r = CubicSpline(thetas, rs, bc_type='periodic')
    #nphis should be increased if the test fails
    nphis = 1000
    phis = np.linspace(0.0, 2.0*np.pi, nphis, endpoint=True)
    ts = np.zeros(nphis)
    ts = integrate(phis, r)
    tmax = ts[-1]
    phi = CubicSpline(ts, phis)
    phi = shift_input(phi, 0.0, tmax)
    test_ts = np.linspace(0.0, tmax/8.0, 1000)
    return np.allclose(phi(test_ts), gudermannian(test_ts/w))

##########################################################################################

def main(n=4, w=1.0, nframes=100, ncycles=6, reverse=False, test=False):
    '''
    Produce images of a shape rolling without slipping on a curve
    n = number of polygon sides
    w = width of polygon
    nframes = number of images to save
    ncycles = number of complete revolutions
    reverse = toggle for reversing the motion
    test = toggle for testing quadrature accuracy
    '''
    #thetas, rs = ellipse()
    thetas, rs = polygon(n, w)
    height = max(rs)
    #r = interpolated function for the polar representation r(theta)
    r = CubicSpline(thetas, rs, bc_type='periodic')

    if test == True:
        close = test_quadrature()
        print('Quadrature result agrees with exact result for square:', close)

    #nphis = number of quadrature points
    nphis = 1000
    phis = np.linspace(0.0, 2.0*np.pi, nphis, endpoint=True)
    ts = integrate(phis, r)
    tmax = ts[-1]

    #phi = interpolated function for the rotation angle phi(t) at time t
    phi = CubicSpline(ts, phis)
    #make phi periodic in time (period tmax)
    phi = shift_input(phi, 0.0, tmax)

    t0 = -ncycles*tmax/2.0
    tf = ncycles*tmax/2.0
    #plot_ts = times at which images will be generated
    #need to make sure we do not redefine ts here because CubicSpline has a reference
    plot_ts = np.linspace(t0, tf, nframes)
    if reverse == True:
        plot_ts = plot_ts[::-1]
    #dense_plot_ts = times used for plotting the curve the shape rolls on
    dense_plot_ts = np.linspace(t0, tf, 100*nframes)
    for j, t in enumerate(plot_ts):
        qs = r(phi(t) + thetas)
        qs = polar_to_cartesian(thetas, qs, np.array([t, 0.0]))
        plt.plot(qs[:, 0], qs[:, 1], lw=3, color='r')
        plt.plot(dense_plot_ts, -r(phi(dense_plot_ts)), lw=3, color='k')
        ax = plt.gca()
        #make sure the square looks square
        ax.set_aspect('equal', adjustable='box')
        plt.xticks([])
        plt.yticks([])
        plt.xlim(t0 + tmax/n, tf - tmax/n)
        plt.ylim(-1.06*height, 1.06*height)
        ax.axis('off')
        fig = plt.gcf()
        fig.set_size_inches(20, 10)
        plt.savefig('figures/%.4i.png' % j, bbox_inches="tight", pad_inches=0, dpi=100)
        plt.clf()

##########################################################################################

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--n', type=int, default=4, help='Number of polygon sides')
  parser.add_argument('--w', type=float, default=1.0, help='Width of polygon')
  parser.add_argument('--nframes', type=int, default=100, help='Number of images to save')
  parser.add_argument('--ncycles', type=int, default=12,
                      help='Number of complete revolutions')
  parser.add_argument('--reverse', action='store_true', default=False,
                      help='Toggle for reversing the motion')
  parser.add_argument('--test', action='store_true', default=False,
                      help='Toggle for testing quadrature accuracy')
  args = parser.parse_args()
  main(args.n, args.w, args.nframes, args.ncycles, args.reverse, args.test)
