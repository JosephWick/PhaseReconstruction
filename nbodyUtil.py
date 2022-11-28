#
# util.py
# Joseph Wick, 2022
#
# static utility functions for P115 project so they don't
# need to be copied between notebooks
#

import numpy as np

class util:

    # nbody()
    # handles math for the nbody simulation when plugged into odeint()
    @staticmethod
    def nbody(Y, t, Ms, isDynamic=True):
        '''
        nbody(...)
            Returns the new state vector for an nbody system. Assumes units such that
            GM = 1, where M is the largest mass in the system

        Parameters
        ----------
            Y : [double]
                Existing state vector. For each body, the state vector will
                have four elements in the following order:
                [x-position, y-position, x-velocity, y-velocity]
            Ms : [double]
                Masses of the bodies. Must satisfy the condition 4*len(Ms) == len(Y)

            isDynamic : [bool]
                ith value is True if that mass moves under gravity. False if that mass
                stays in place and only affects other masses. Defaults to all masses
                being dynamic

        Returns
        -------
            newY : [double]
                Updated state vector following the same convention as parameter Y

        '''
        N = len(Ms)

        # double check that sizes are in convention
        ndof = 4
        if (N*ndof != len(Y)):
            print('nbody: Y must have 4 fields for each body')
            return

        # check if all masses are dynamic
        dyn = []
        if not isDynamic:
            dyn = isDynamic
        else:
            for i in range(N): dyn.append(True)

        newY = np.zeros(len(Y))

        # for each body, update state vector
        # i will iterate over receivers, j over sources
        for i in range(N):
            if not dyn[i]: continue
            # some indices for easier coding
            x  = (ndof*i)
            y  = (ndof*i) +1
            vx = (ndof*i) +2
            vy = (ndof*i) +3

            # for use in acceleration equations
            recx = Y[x]
            recy = Y[y]

            # update position and velocity based on values from the last time step
            newY[x] = Y[vx]
            newY[y] = Y[vy]

            # calculate new accelerations
            # forces are from each other body, so accelerations are as well
            for j in range(N):
                # no gravity on self
                if i==j:
                    continue

                # define quantities
                srcM = Ms[j]
                srcx = Y[ndof*j]
                srcy = Y[ndof*j+1]
                r = ((recx-srcx)**2 + (recy-srcy)**2 )**(3/2)

                # calculations
                newY[vx] += -srcM*(recx-srcx)/r
                newY[vy] += -srcM*(recy-srcy)/r

        return newY

    # nbody()
    # handles math for the nbody simulation when plugged into odeint()
    # first two bodies only interact with eachother (binary)
    # others only feel forces from binary
    # for the purpose of simulating multiple starting positions for a planet at once
    @staticmethod
    def nbodyMult(Y, t, Ms, isDynamic=True):
        '''
        nbody(...)
            Returns the new state vector for an nbody system. Assumes units such that
            GM = 1, where M is the largest mass in the system.

            First two bodies in Y only interact with eachother. The other bodies only
            feel forces from the binary.

        Parameters
        ----------
            Y : [double]
                Existing state vector. For each body, the state vector will
                have four elements in the following order:
                [x-position, y-position, x-velocity, y-velocity]
            Ms : [double]
                Masses of the bodies. Must satisfy the condition 4*len(Ms) == len(Y)

            isDynamic : [bool]
                ith value is True if that mass moves under gravity. False if that mass
                stays in place and only affects other masses. Defaults to all masses
                being dynamic

        Returns
        -------
            newY : [double]
                Updated state vector following the same convention as parameter Y

        '''
        N = len(Ms)

        # double check that sizes are in convention
        ndof = 4
        if (N*ndof != len(Y)):
            print('nbodyMult: Y must have 4 fields for each body')
            return

        # check at least a binary
        if N < 2:
            print('nbodyMult: must be at least 2 bodies')
            return

        # check if all masses are dynamic
        dyn = []
        if not isDynamic:
            dyn = isDynamic
        else:
            for i in range(N): dyn.append(True)

        newY = np.zeros(len(Y))

        # handle interactions of binary
        for i in range(2):
            j=1
            if i==1: j=0
            x  = (ndof*i)
            y  = (ndof*i) +1
            vx = (ndof*i) +2
            vy = (ndof*i) +3

            newY[x] = Y[vx]
            newY[y] = Y[vy]

            # for use in acceleration equations
            recx = Y[x]
            recy = Y[y]

            srcM = Ms[j]
            srcx = Y[ndof*j]
            srcy = Y[ndof*j+1]
            r = ((recx-srcx)**2 + (recy-srcy)**2 )**(3/2)

            newY[vx] += -srcM*(recx-srcx)/r
            newY[vy] += -srcM*(recy-srcy)/r

        # handle interactions of planets
        # i is receiver, j is source
        for i in range(2,N):
            if not dyn[i]: continue
            # some indices for easier coding
            x  = (ndof*i)
            y  = (ndof*i) +1
            vx = (ndof*i) +2
            vy = (ndof*i) +3

            # for use in acceleration equations
            recx = Y[x]
            recy = Y[y]

            # update position and velocity based on values from the last time step
            newY[x] = Y[vx]
            newY[y] = Y[vy]

            # calculate new accelerations
            # forces are from each other body, so accelerations are as well
            for j in range(2):
                # no gravity on self
                if i==j:
                    continue

                # define quantities
                srcM = Ms[j]
                srcx = Y[ndof*j]
                srcy = Y[ndof*j+1]
                r = ((recx-srcx)**2 + (recy-srcy)**2 )**(3/2)

                # calculations
                newY[vx] += -srcM*(recx-srcx)/r
                newY[vy] += -srcM*(recy-srcy)/r

        return newY

    @staticmethod
    def getAxisPlanets(x,y,vx,vy):
        '''
        getAxisPlanets(...)
            Given x,y, vx,vy for a planet that starts with y=0, returns state
            vector for four planets with same r and v

        Parameters
        ----------
        x : double
            x position of planet

        y : double
            y position of planet

        vx : double
            x velocity of planet

        vy : double
            y velocity of planet

        Returns
        -------
        Y : [double]
            state vector for four planets. For each planet there are four elements in
            the order [xposition, yposition, xvelocity, yvelocity]

        '''

        Y = [x,y, vx, vy]
        Y = Y + [y,x, -vy, vx]
        Y = Y + [-x,y, -vx,-vy]
        Y = Y + [y,-x, vy, -vx]

        return Y

    # getBin()
    # returns intitial conditions for a specific binary
    @staticmethod
    def getBin():
        '''
        getBin()
            returns a list of the masses of the binary and a list of the initial conditions

        Returns
        -------
        Ms : List
            List of the masses of the binary, normalized such that G*M = 1

        Y : List
            List of the initial conditions of each mass of the binary

        '''
        # masses
        M = 1
        m = 0.67

        # used to counter movement of the entire system in one direction
        voff = -0.5817365
        xoff = -0.3

        # larger mass
        xM  = 0. + xoff
        yM  = 0.
        vxM = 0.0
        vyM = 0.0 + voff
        # smaller mass
        xm  = 0.75 + xoff
        ym  = 0.
        vxm = 0.0
        vym = 1.45 + voff

        # assemble
        Ms = [M,m]
        Y = [xM,yM,vxM,vyM, xm,ym,vxm,vym]

        return Ms, Y

    # HW()
    # Holman & Wiegart criterion
    @staticmethod
    def HW(e,mu, binsep):
        '''
        HW(...)
            calculates the upper and lower limit of the Holman & Wiegart criterion for stability of orbits around a binary

        Parameters
        ----------
        e : double
            Eccentricity of the binary

        mu : double
            Mass ratio of the binary. Defined as mu = ma/(ma+mb) where ma is the larger mass.

        Returns
        -------
        HWp : Double
            upper limit semi-major axis of a stable orbit in units of binary separation

        HWm : Double
            lower limit semi-major axis of a stable orbit in units of binary separation
        '''
        HWp = 1.64 + 5.15*e + -2.11*(e**2) + 4.21*mu + -4.10*e*mu + -4.98*(mu**2) + 4.97*(e**2)*(mu**2)

        HWm = 1.56 + 5.05*e + -2.33*(e**2) + 4.03*mu + -4.44*e*mu + -5.20*(mu**2) + 4.25*(e**2)*(mu**2)

        return HWp*binsep, HWm*binsep

    # getEscapeV()
    # returns escape velocity as a function of distance and central mass
    @staticmethod
    def getEscapeV(M,r):
        '''
        only works if G*M is defined to be 1
        '''
        return np.sqrt(2*M / r)


    # getUCV()
    # return uniform circular motion velocity
    @staticmethod
    def getUCV(M,r):
        '''
        getUCV(...)
            Calculates the velocity for uniform circular motion around a body of mass M
            at radius r.

        Parameters
        ----------
            M : double
                Mass of central mass

            r : double
                radius at which body of concern moves around

        Returns
        -------
            UCV : double
                Velocity for uniform circular motion
        '''

        return np.sqrt(M/r)

    # getUCP()
    # returns period of uniform circular motion
    @staticmethod
    def getUCP(r,v):
        return 2*np.pi*r / v

    # getComponents()
    # returns x and y components given a magnitude and a theta
    @staticmethod
    def getComponents(mag, theta):
        '''
        getComponents(...)
            returns x and y components given a magnitude and a theta

        Parameters
        ----------
        mag : double
            magnitude

        theta : double
            direction

        Returns
        -------
        x : double
            xhat component

        y : double
            yhat component
        '''
        return mag*cos(theta), mag*sin(theta)

    # isStable()
    # returns true if an integrated state vector is stable, false otherwise
    @staticmethod
    def isStable(Y):
        '''
        isStable(...)
            boolean check for if a planet is stable under two conditions:
                1. velocity never is greater than escape velocity
                2. planet never passes interior to the binary orbit
            Currently this assumes two stars and one planet in the order M,m,p

        Parameters
        ----------
        Y : [double]
            state vector of the system over the duration of the simulation

        Returns
        -------
        stability : bool
            True if the planet is stable. False otherwise
        '''

        # check v and vmax
        r = np.sqrt(Y[0,8]**2 + Y[0,9]**2)
        vmax = util.getEscapeV(1.67,r)
        v = np.sqrt(Y[:,10]**2 + Y[:,11]**2)
        if v.max() > vmax:
            #print('v')
            return False

        # check if planet goes in binary separation
        # if at any point r<0.5 from center of binary
        bcenx = Y[:,4].mean()
        bceny = Y[:,5].mean()
        r = np.sqrt((Y[:,8]-bcenx)**2 + (Y[:,9]-bceny)**2)
        if r.min() < 0.5:
            #print('r')
            return False

        return True
