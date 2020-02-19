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
        self.correction = 0  # for double integral factor to find deflection
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

    def y_tag_tag(self, position):  # E*I*y''=M
        return (self.moment(position) / (self.e * self.i))

    def y_tag(self, position):
        res = scipy.integrate.quad(self.y_tag_tag, 0, position)+self.correction
        return(res[0])

    def y(self, position):
        res = scipy.integrate.quad(self.y_tag, 0, position)
        return(res[0])

    def correct_factor(self):
        #find defelection at right support
        test = self.y(self.right)
        if test > 1:
            count = 0
            while count < 100 and test > 1:
                self.correction -= 0.2
                test = self.y(self.right)
                count+=1
            if count < 100:
                print('Converged! y is', test,'at right support')
                exit()
            print('not converged while test>1', test)

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
        input('press enter...')


f1 = linear_load(0, -270, 0, 6)
f2 = force(-600, 3)
f3 = force(-900, 4)
cl1 = constant_load(-60, 0, 18)
m = moment(4800, 9)
a = simple_support_2(5, [f2, f3], e=1, i=1)
print('reactions:', a.r1.magnitude(), a.r2.magnitude())
print('moment at', a.length, a.moment(a.length))
print('deflection at 2.5', a.y(2.5))
a.plot_moment()
