import math
import scipy.integrate

import matplotlib.pyplot as plt
plt.ion()


class force:
    type = 'force'

    def __init__(self, magnitude=0, position=0):
        self.mag = magnitude
        self.pos = position
        # where the load is acting, to be the same definition as linear loads
        self.start = position

    def position(self):
        return (self.pos)

    def magnitude(self):
        return(self.mag)


class moment(force):
    type = 'moment'


class constant_load:
    type = 'constant_load'

    def __init__(self, load=0, start_pos=0, end_pos=0):
        self.load = load  # magnitude units is N/m
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.start = start_pos

    def position(self):
        return ((self.end_pos + self.start_pos) / 2)

    def magnitude(self):
        return(self.load*(self.end_pos-self.start_pos))


class linear_load:
    type = 'linear_load'

    def __init__(self, start_mag=0, end_mag=0, start_pos=0, end_pos=0):
        self.start_mag = start_mag
        self.end_mag = end_mag
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.start = start_pos
        if (end_pos - start_pos) < 1e-6:
            self.slope = 0
        else:
            self.slope = math.atan((end_mag-start_mag)/(end_pos-start_pos))

    def magnitude(self):
        return((self.start_mag+self.end_mag)*(self.end_pos-self.start_pos)/2)

    def position(self):
        nominator = (self.end_pos - self.start_pos) * \
            (self.start_mag + 2 * self.end_mag)
        denominator = 3 * (self.start_mag + self.end_mag)
        if denominator != 0:
            return (self.start_pos+(nominator / denominator))
        else:
            return (0)


class simple_support_2:  # a simple beam with simple support on both of its ends
    def __init__(self, length=1, loads=[], left=0, right=0, e=70e9, i=1):
        self.length = length
        self.loads = loads
        self.e = e
        self.i = i
        self.correction = 1  # for double integral factor to find deflection
        self.correction2 = 0
        for load in self.loads:
            if load.position() > self.length:
                print('Error - one of the loads exceeds beam length!')
                exit()
        self.left = left
        if right == 0:
            self.right = self.length
        else:
            self.right = right
        self.calculate_reactions()

    def calculate_reactions(self):  # should only be run once in initialization!
        sum_moments = 0
        sum_magnitude = 0
        for load in self.loads:
            if load.type != 'moment':
                if load.position() < self.left:
                    sum_moments += load.magnitude() * (self.left-load.position())
                else:
                    sum_moments -= load.magnitude()*(load.position()-self.left)
                sum_magnitude += load.magnitude()
            else:
                sum_moments -= load.magnitude()
        self.r2 = force(sum_moments/(self.right-self.left), self.right)
        self.r1 = force(-1 * (self.r2.magnitude() + sum_magnitude), self.left)
        self.loads += [self.r1, self.r2]

    # calculate shear force at position from left support. position must be between 0 and length

    def __shear_moment(self, position):
        if position > self.length:
            print('Error - position exceeds beam length!')
            exit()
        sum_shears = 0
        sum_moments = 0
        for load in self.loads:
            if load.start <= position:
                if load.type == 'linear_load':
                    if position < load.end_pos:
                        end_mag = load.start_mag + \
                            (position-load.start)*math.tan(load.slope)
                        load = linear_load(
                            load.start_mag, end_mag, load.start, position)
                if load.type == 'constant_load':
                    if position < load.end_pos:
                        load = constant_load(
                            load.load, load.start_pos, position)
                if load.type != 'moment':
                    sum_shears += load.magnitude()
                    sum_moments += load.magnitude() * load.position()
                else:
                    sum_moments += load.magnitude()
        return ([sum_shears, sum_moments])

    def shear(self, position):
        temp = self.__shear_moment(position)
        return (temp[0])

    def moment(self, position):
        temp = self.__shear_moment(position)
        return (-1 * temp[1] + position * temp[0])

    def moment_area(self, a, b):  # compute moment area between b>=x>=a
        if (a > b or a < 0 or a > self.length or b < 0 or b > self.length):
            print('error - index error in moment_area')
            exit()
        return (scipy.integrate.quad(self.moment, a, b)[0])

    def moment_da_x(self, position):  # function needed for moment_x_cg integration
        return(self.moment(position)*position)

    def moment_x_cg(self, a, b):  # compute the x_cg of the moment from a over span a to b
        f = scipy.integrate.quad(self.moment_da_x, a, b)[0]
        return (f / self.moment_area(a, b))

    def y(self, position):  # compute defelction using area-moment theorem
        # a is left support, b is position to find, c is right support
        # refer to https://www.mathalino.com/reviewer/mechanics-and-strength-of-materials/deflection-in-simply-supported-beams-area-moment-method

        if position > self.left and position < self.right:
            # find xc, which is the distance to cg of the moment area, from support c
            xc = self.moment_x_cg(self.left, self.right)
            xc = (self.right - self.left) - xc
            # now compute t_ca
            t_ca = self.moment_area(
                self.left, self.right) * xc / (self.e * self.i)
            # find xb, which is the distance to cg of the moment area, from position to support a
            xb = self.moment_x_cg(self.left, position)
            xb = (position - self.left) - xb
            # now compute t_ba
            t_ba = self.moment_area(self.left, position) * \
                                    xb / (self.e * self.i)
            return(-1*(position*t_ca/(self.right-self.left)-t_ba))

        if position == self.left or position == self.right:
            return (0)

        if position > self.right:
            xa = self.moment_x_cg(self.left, self.right)
            xa = xa - self.left
            t_ab = self.moment_area(
                self.left, self.right) * xa / (self.e * self.i)
            print('t_ab is', t_ab)
            y_position = t_ab * (position - self.right) / \
                                 (self.right - self.left)
            xc = self.moment_x_cg(self.right, position)
            xc = position - xc
            print('xc is', xc)
            print('moment_area bc is', self.moment_area(self.right, position))
            t_cb = self.moment_area(
                self.right, position) * xc / (self.e * self.i)
            print('t_cb is', t_cb)
            return (y_position + t_cb)

        if position < self.left:
            xa = self.moment_x_cg(position, self.left)
            xa = xa - position
            t_ab = self.moment_area(position, self.left) * \
                                    xa / (self.e * self.i)
            xc = self.moment_x_cg(self.left, self.right)
            xc = self.right - xc
            t_cb = self.moment_area(
                self.left, self.right) * xc / (self.e * self.i)
            y_position = t_cb * (self.left - position) / \
                                 (self.right - self.left)
            return(y_position+t_ab)

    def plot_moment(self):
        data_x = []
        data_y = []
        x = 0
        delta = 0.1
        while x <= self.length:
            data_x += [x]
            data_y += [self.moment(x)]
            x += delta
        plt.plot(data_x, data_y)
        return

    def plot_deflection(self):
        data_x = []
        data_y = []
        x = 0
        delta = 0.1
        while x <= self.length:
            data_x += [x]
            data_y += [self.y(x)]
            x += delta
        plt.plot(data_x, data_y)
        return

    


f1 = linear_load(0, -900, 1, 4)
m1=moment(600,0)
f2 = force(-60, 13)
f3 = force(400, 7)
cl1 = constant_load(-80, 0, 4)
m = moment(4800, 9)
a = simple_support_2(4, [f1,m1], e=1, i=1,left=1)
print('reactions:', a.r1.magnitude(), a.r2.magnitude())
print('moment at', a.length, a.moment(a.length))
print('deflection at 0:', a.y(4))

a.plot_moment()
a.plot_deflection()
input('press enter...')
