import math
from schauderhat2mv_helper import searchMinkowski, simpDenom, checkUni, interpolate
from config import Symbols
from sigmaConstruct import sigmaConstruct, sigmaConstruct_agu


def subdivide(cords): #subdivide the simplicial convex [n1/d1, n2/d2] until it is unimodular
    cords = simpDenom(cords)
    if all(checkUni(cords)): return cords
    Unicheck = [False] #keep track of if all simplices are unimodular
    while not all(Unicheck):
        indices = [i for i, x in enumerate(Unicheck) if x == False] #find all that are not unimodular
        addvertices = []
        for idx in indices:
            (u,v) = searchMinkowski(cords[idx], cords[idx+1])
            addvertices.append((u, v))

        #now insert addvertices back to cords
        assert (len(addvertices) == len(indices))
        acc = 1
        for i in range(len(addvertices)):
            cords.insert(indices[i] + acc, addvertices[i])
            acc += 1
        #again check if unimodular
        cords = simpDenom(cords)
        Unicheck = checkUni(cords)
    return cords

cords = [(0,1), (1,4), (2, 4), (3,4), (1,1)]


cords = simpDenom(cords)
x=2


class schauderhat2mv:
    def __init__(self, xcords, ycords):  # input: xcoordinates, ycoordinates of the breakpoints
        self.xcords = xcords
        self.ycords = ycords

    def subdivide(self):
        cords = self.xcords.copy()
        cords = simpDenom(cords)
        Unicheck = checkUni(cords) # keep track of if all simplices are unimodular
        # if all(Unicheck): return cords
        while not all(Unicheck):
            indices = [i for i, x in enumerate(Unicheck) if x == False]  # find all that are not unimodular
            addvertices = []
            for idx in indices:
                (u, v) = searchMinkowski(cords[idx], cords[idx + 1])
                addvertices.append((u, v))

            # now insert addvertices back to cords
            assert (len(addvertices) == len(indices))
            acc = 1
            for i in range(len(addvertices)):
                cords.insert(indices[i] + acc, addvertices[i])
                acc += 1
            # again check if unimodular
            cords = simpDenom(cords)
            Unicheck = checkUni(cords)
        self.univertices = cords

    def buildhats(self):
        self.hats = [None] * (len(self.univertices) - 2) #list of hats, each hat is [(w1, b1), (w2, b2)]. The left piece is y = w1x+b1, the right piece is y = w2x+b2
        for i in range(len(self.hats)):
            n1, d1 = self.univertices[i]
            n2, d2 = self.univertices[i+1] #the vertex in the middle
            n3, d3 = self.univertices[i+2]
            wl = d1 // (d1*n2 - d2*n1)
            bl = n1 // (d2*n1 - d1*n2)
            wr = d3 // (d3*n2 - d2*n3)
            br = n3 // (d2*n3 - d3*n2)
            self.hats[i] = [(wl, bl), (wr, br)]
    def convert(self, method):
        #first convert to normal form
        cordscopy = self.xcords.copy()
        cordscopy = simpDenom(cordscopy)
        #calculate the magnitude of each hat
        mags = [0] * len(self.univertices)
        indices = [0] * len(cordscopy)
        for i in range(1, len(cordscopy)):
            indices[i] = self.univertices.index(cordscopy[i])
        for i in range(len(indices)):
            idx = indices[i]
            mags[idx] = self.ycords[i] * self.univertices[idx][1]
        for i in range(len(self.xcords)-1):
            if indices[i+1] - indices[i] == 1:
                continue
            for j in range(indices[i]+1, indices[i+1]):
                mags[j] = interpolate(xl=self.xcords[i], yl=self.ycords[i], xr=self.xcords[i+1], yr=self.ycords[i+1], mid = self.univertices[j])
        self.hatmags = mags[1:-1]

        #find the mv term of each hat function

        self.mvhat = [''] * len(self.hats)
        for i in range(len(self.hats)):
            left = sigmaConstruct_agu(self.hats[i][0][0], self.hats[i][0][1])
            right = sigmaConstruct_agu(self.hats[i][1][0], self.hats[i][1][1])
            self.mvhat[i] = f"({left}) {Symbols.WEDGE} ({right})"



    def compose(self):
        self.mvhatcomp = [''] * len(self.hats)
        for i in range(len(self.hats)):
            if self.hatmags[i] == 0:
                continue
            self.mvhatcomp[i] = Symbols.OR.join([f"h{i}"]*self.hatmags[i])
        self.overallmv = Symbols.OR.join(self.mvhatcomp)
        for i in range(len(self.hats)):
            if self.hatmags[i] == 0:
                continue
            self.overallmv = self.overallmv.replace(f"h{i}", self.mvhat[i])





    def extract(self, method):
        self.subdivide() #subdivide until unimodular
        self.buildhats() #construct a schauder hat on each vertix
        self.convert(method) #convert to mv terms
        self.compose() #compose to get an overall mv term

