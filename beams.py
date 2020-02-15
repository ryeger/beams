import math


class force:
    type = 'force'

    def __init__(self, magnitude=0, position=0):
        self.magnitude = magnitude
        self.position = position
        # where the load is acting, to be the same definition as linear loads
        self.start = position


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
        return ((self.end_pos - self.start_pos) / 2)

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
        self.slope = (end_mag-start_mag)/(end_pos-start_pos)

    def magnitude(self):
        return((self.start_mag+self.end_mag)*(self.end_pos-self.start_pos)/2)

    def position(self):
        nominator = (self.end_pos - self.start_pos) * \
            (self.start_mag + 2 * self.end_mag)
        denominator = 3 * (self.start_mag + self.end_mag)
        if denominator != 0:
            return (nominator / denominator)
        else:
            return (0)


class simple_support_2:  # a simple beam with simple support on both of its ends
    def __init__(self, length=1, loads=[], e=70e9, i=1):
        self.length = length
        self.loads = loads
        for load in self.loads:
            if load.position > self.length:
                print('Error - one of the loads exceeds beam length!')
                exit()
        self.calculate_reactions()

    def calculate_reactions(self):
        sum_moments = 0
        sum_magnitude = 0
        for load in self.loads:
            if load.type != 'moment':
                sum_moments += load.magnitude * load.position
                sum_magnitude += load.magnitude
            else:
                sum_moments += load.magnitude
        self.r2 = (-1 / self.length) * sum_moments
        self.r1 = -1*(self.r2+sum_magnitude)

    def shear(self, position):  # calculate shear force at position from left support. position must be between 0 and length
        if position > self.length:
            print('Error - position exceeds beam length!')
            exit()
        sum_shear = self.r1
        for load in self.loads:
            if load.start <= position:
                if load.type == 'linear_load':
                    end_mag = load.start_mag + \
                        (position-load.start)*math.tan(load.slope)
                    load = linear_load(load.start_mag,end_mag,load.start,position)
                if load.type == 'constant_load':
                    load = constant_load(load.load, load.start_pos, position)
                if load.type != 'moment':
                    sum_shear += load.magnitude
        return(-1*sum_shear)

f = force(40, 10)
m = moment(50, 12)
a = simple_support_2(20, [f])
print('reactions:', a.r1, a.r2)
print('shear at position 11',a.shear(11))
