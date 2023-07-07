from cmath import pi
import math
from pybbat.bmath import bmath
import numpy as np


class bTools:

    ###
    # accmath
    ###
    def qguass_1t(func, a, b, t, t1):
        """Performs numerical integeration using the Gauss-Legendre quadrature formula with 10 points
            Evaluates the integral of the provided 'func' over the interval [a,b] by first mapping the interval [-1,1]
            The transformation used is x = xm + xr*t where xm= 9b+a)/2 and xr = (b-a)/2
            Approximates the integral using the sum of weighted function evals at the Gauss-Legendre quadrature
            points in arrays x and w

        Args:
            func (function): function to integrate
            a (float): minimum of the interval to integrate over
            b (float): maximum of the interval to integrate over
            t (float): inital time (for the integration)
            t1 (float): final time (for the integration)
        """
        npoint = 10
        xr = 0.5 * (b - a)
        xm = 0.5 * (b + a)
        s = 0
        x = [
            0.0,
            0.076526521133497,
            0.227785851141645,
            0.373706088715419,
            0.510867001950827,
            0.636053680726515,
            0.746331906460150,
            0.839116971822218,
            0.912234428251325,
            0.963971927277913,
            0.993128599185094,
        ]
        w = [
            0.0,
            0.152753387130725,
            0.149172986472603,
            0.142096109318382,
            0.131688638449176,
            0.118194531961518,
            0.101930119817240,
            0.083276741576704,
            0.062672048334109,
            0.040601429800386,
            0.017614007139152,
        ]

        for i in range(0, npoint):
            dx = xr * x[i]
            s += w[i] * (func(xm + dx, t, t1) + func(xm - dx, t, t1))
        return s * xr

    def rtsafe_l_t(funcd, a, b, alpha, phi_s):
        """Secant Method implementation
            finds the root of the function in a given interval

        Args:
            funcd (function): function to find the root of
            a (float): minimum of the interval
            b (float): maximum of the interval
            alpha (float): alpha value for the function
            phi_s (float): phi value for the function

        Returns:
            float: rts when the iterations stop and inf if there is an error/never converges
        """
        # domain (a,b) f(a)<0,f(b)>0
        MAXIT = 200
        epsilon = 1e-14
        xl = a
        xh = b
        rts = 0.5 * (a + b)
        dxold = abs(b - a)
        dx = dxold
        f, df = funcd(rts, alpha, phi_s)
        for j in range(0, MAXIT):
            if (((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0) or (
                abs(2.0 * f) > abs(dxold * df)
            ):
                dxold = dx
                dx = 0.5 * (xh - xl)
                rts = xl + dx
                if xl == rts:
                    return rts
            else:
                dxold = dx
                dx = f / df
                temp = rts
                rts -= dx
                if temp == rts:
                    return rts
            if abs(dx) < epsilon:
                return rts
            f, df = funcd(rts, alpha, phi_s)
            if f < 0.0:
                x1 = rts
            else:
                xh = rts
        return 1
        # return float("inf")

    def qchebyshev_t(func, a, b, t):
        """Computes the chebyshev using 20 points and given coefficients

        Args:
            func (function): python function representing the numerical equation
            a (float): float input to func
            b (float): float input to func
            t (float): float input to func

        Returns:
            float : value from chebyshev approximation
        """
        npoint = 20
        x = [
            0.0,
            0.996917333733127,
            0.972369920397676,
            0.923879532511286,
            0.852640164354092,
            0.760405965600030,
            0.649448048330183,
            0.522498564715948,
            0.382683432365089,
            0.233445363855905,
            0.078459095727845,
            -0.078459095727844,
            -0.233445363855905,
            -0.382683432365089,
            -0.522498564715948,
            -0.649448048330183,
            -0.760405965600030,
            -0.852640164354092,
            -0.923879532511286,
            -0.972369920397676,
            -0.996917333733127,
        ]

        xm = 0.5 * (b + a)
        xr = 0.5 * (b - a)
        s = 0
        for i in range(0, npoint):
            yj = xm + xr * x[i]
            s += func(yj, a, b, t)
        return s * bmath.pi / npoint

    def qchebyshev_7t(func, vrf, vn, n, theta, t, a, b):
        """Computes the chebyshev using 20 points and given coefficients

        Args:
            func (function): python function representing the numerical equation
            vrf (float): RF voltage/turn (kV)
            vn (float): Vn value
            n (float): number of bunches
            theta (float): theta in degrees
            t (float): float input to func
            a (float): float input to func
            b (float): float input to func

        Returns:
            float : value from chebyshev approximation
        """
        npoint = 20
        x = [
            0.0,
            0.996917333733127,
            0.972369920397676,
            0.923879532511286,
            0.852640164354092,
            0.760405965600030,
            0.649448048330183,
            0.522498564715948,
            0.382683432365089,
            0.233445363855905,
            0.078459095727845,
            -0.078459095727844,
            -0.233445363855905,
            -0.382683432365089,
            -0.522498564715948,
            -0.649448048330183,
            -0.760405965600030,
            -0.852640164354092,
            -0.923879532511286,
            -0.972369920397676,
            -0.996917333733127,
        ]

        xm = 0.5 * (b + a)
        xr = 0.5 * (b - a)
        s = 0
        for j in range(0, npoint):
            yj = xm + xr * x[j]
            s += func(yj, vrf, vn, n, theta, t, a, b)

        return s * bmath.pi / npoint

    def qgauss_7t(func, vrf, vn, n, theta, phis, phi1, phi2):
        """Computes the gaussian using 10 points and given coefficients

        Args:
            func (function): python function representing the numerical equation
            vrf (float): RF voltage/turn (kV)
            vn (float): Vn value
            n (float): number of bunches
            theta (float): theta in degrees
            t (float): float input to func
            a (float): float input to func
            b (float): float input to func

        Returns:
            float : value from chebyshev approximation
        """
        npoint = 10
        x = [
            0.0,
            0.076526521133497,
            0.227785851141645,
            0.373706088715419,
            0.510867001950827,
            0.636053680726515,
            0.746331906460150,
            0.839116971822218,
            0.912234428251325,
            0.963971927277913,
            0.993128599185094,
        ]
        w = [
            0.0,
            0.152753387130725,
            0.149172986472603,
            0.142096109318382,
            0.131688638449176,
            0.118194531961518,
            0.101930119817240,
            0.083276741576704,
            0.062672048334109,
            0.040601429800386,
            0.017614007139152,
        ]

        xm = 0.5 * (phi2 + phi1)
        xr = 0.5 * (phi2 - phi1)
        s = 0
        for j in range(0, npoint):
            dx = xr * x[j]
            s += w[j] * func(xm + dx, vrf, vn, n, theta, phis, phi1, phi2) + func(
                xm - dx, vrf, vn, n, theta, phis, phi1, phi2
            )

        return s * xr

    ###
    # bktbun.c
    ###
    def sign(A):
        """Returns the sign of float input A

        Args:
            A (float): value to verify sign of

        Returns:
            int : positive 1 if A is greater than 0 or -1 if A is less than 0
        """
        return 1 if A >= 0 else -1

    def cvector(nl, nh):
        """function creates a dynamic array of characters

        Args:
            nl (int): lower value of the array
            nh (int): upper value of the array

        Returns:
            array: array of characters
        """
        v = bytearray(nh - nl + 1)
        return v[nl - 1 :]

    def ivector(nl, nh):
        """
        Allocate and return an integer vector with indices ranging from nl to nh.

        Args:
            nl (int): The lower bound index of the vector.
            nh (int): The upper bound index of the vector.

        Returns:
            numpy.ndarray: An integer vector with indices ranging from nl to nh.
        """
        v = np.zeros(nh - nl + 1, dtype=int)
        return v - nl

    def dvector(nl, nh):
        """
        Allocates a double vector with subscript range v[nl..nh].

        Parameters:
        nl (int): The lower subscript of the vector.
        nh (int): The upper subscript of the vector.

        Returns:
        list: A list representing the allocated double vector with subscript range v[nl..nh].
        """
        v = (nh - nl + 1) * [0.0]
        return v

    def proper_phi(phis, phi):
        """
        Convert a RF generic angle into its proper position according to phis.

        Args:
            phis (float): A float representing the value of phis.
            phi (float): A float representing the value of phi.

        Returns:
            A float representing the proper position of phi.
        """
        tmphis = phis

        if 0 <= tmphis < bmath.pi / 2:
            return phi
        elif bmath.pi / 2 <= tmphis < bmath.pi:
            return bmath.pi - phi
        elif bmath.pi <= tmphis < bmath.pi + bmath.pi / 2:
            return bmath.pi + phi
        elif bmath.pi + bmath.pi / 2 <= tmphis <= bmath.pi + bmath.pi:
            return 2 * bmath.pi + phi

    def phi_s(etas, C, rho, Bdot, Vrf):
        """
        Calculates the RF phase angle phi_s according to the given parameters.

        Parameters:
        etas (float): the synchronous phase angle
        C (float): the machine constant
        rho (float): the particle's mass-to-charge ratio
        Bdot (float): the time derivative of the magnetic field
        Vrf (float): the RF voltage

        Returns:
        float: the calculated RF phase angle phi_s
        """
        tmp = C * rho * Bdot / Vrf
        # print("phi_s")
        # print(C)
        # print(rho)
        # print(Bdot)
        # print(Vrf)
        # print(tmp)
        # print("end")

        if abs(tmp) > 1:
            return float("inf")

        if Vrf > 0.0:
            # print("if")
            return math.asin(tmp) if etas < 0 else bmath.pi - math.asin(tmp)
        else:
            # print("else")
            return (
                2 * bmath.pi - math.asin(tmp) if etas < 0 else bmath.pi + math.asin(tmp)
            )

    def phi_2_bkt(phis):
        """
        Convert the angle phis to the corresponding angle in the second bucket.

        Args:
        phis: A float representing the angle phis.

        Returns:
        A float representing the angle in the second bucket.
        """
        if phis >= 0 and phis < bmath.pi_2:
            return bmath.pi - phis
        if phis >= bmath.pi_2 and phis <= bmath.pi:
            return bmath.pi - phis
        if phis >= bmath.pi and phis < bmath.pi + bmath.pi_2:
            return 3 * bmath.pi - phis
        if phis >= bmath.pi + bmath.pi_2 and phis <= bmath.pi + bmath.pi:
            return 3 * bmath.pi - phis

    def generic_phis(phis: float) -> float:
        """
        Convert a generic angle phis to its corresponding angle within the range [0, 2*pi).

        Args:
        phis: A float representing the generic angle phis.

        Returns:
        A float representing the corresponding angle within the range [0, 2*pi).
        """
        tmphis = phis

        if 0 <= tmphis < bmath.pi / 2:
            return tmphis
        elif bmath.pi / 2 <= tmphis < bmath.pi:
            return bmath.pi - tmphis
        elif bmath.pi <= tmphis < bmath.pi + bmath.pi / 2:
            return tmphis - bmath.pi
        elif bmath.pi + bmath.pi / 2 <= tmphis <= bmath.pi + bmath.pi:
            return 2 * bmath.pi - tmphis

    def generic_phi2(phis, phi2):
        """Computes the phi2 value

        Args:
            phis (float): current value of phis
            phi2 (float): current value of phis2

        Returns:
            float: new value of phi2
        """
        if phis >= 0 and phis < bmath.pi_2:
            return phi2

        if phis >= bmath.pi_2 and phis <= bmath.pi:
            return bmath.pi - phi2

        if phis >= bmath.pi and phis < (bmath.pi + bmath.pi_2):
            return phi2 - bmath.pi

        if phis >= bmath.pi + bmath.pi_2 and phis <= (bmath.pi + bmath.pi):
            return 2 * bmath.pi - phi2

    def bktrtsafe(funcd, phis):
        """
        Find the root of the function funcd using the Brent's method with safe-guarded interval.

        Args:
        funcd: A function that takes a double x and phis, and returns a tuple of double values representing f and df.
        phis: A float representing the value of phis.

        Returns:
        A float representing the estimated root of the function.
        """

        MAXIT = 100
        epsilon = 1e-8

        if phis <= epsilon:
            return -bmath.pi
        if bmath.pi_2 - phis <= epsilon:
            return -bmath.pi_2

        x1 = -bmath.pi - phis
        x2 = phis

        xl = x1
        xh = x2

        rts = 0.5 * (x1 + x2)
        dxold = abs(x2 - x1)
        dx = dxold
        f, df = funcd(rts, phis)
        for j in range(1, MAXIT + 1):
            if (((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0) or (
                abs(2.0 * f) > abs(dxold * df)
            ):
                dxold = dx
                dx = 0.5 * (xh - xl)
                rts = xl + dx
                if xl == rts:
                    return rts
            else:
                dxold = dx
                dx = f / df
                temp = rts
                rts -= dx
                if temp == rts:
                    return rts
            if abs(dx) < epsilon:
                return rts
            f, df = funcd(rts, phis)
            if f < 0.0:
                xl = rts
            else:
                xh = rts

        return float("inf")

    def fnphi1(x, phis):
        """computes the function and derivative at x and phis

        Args:
            x (float): value to integrate over
            phis (float): synchronous phase angle

        Returns:
            float, float: the value at x and phis and the derivative at x and phis
        """
        fn = math.cos(x) + math.cos(phis) + (x - (math.pi - phis)) * math.sin(phis)
        df = math.sin(phis) - math.sin(x)
        return fn, df

    def phi_1_bkt(phis):
        """
        Calculate the value of phi_1_bkt for the given phis.

        Args:
        phis: A float representing the value of phis.

        Returns:
        A float representing the calculated value of phi_1_bkt.
        """

        tmphis = phis
        out = bTools.bktrtsafe(bTools.fnphi1, phis)
        if math.isinf(out):
            return float("inf")

        if tmphis >= 0 and tmphis < math.pi / 2:
            phis = tmphis
            out = bTools.bktrtsafe(bTools.fnphi1, phis)
            return out

        if tmphis >= math.pi / 2 and tmphis <= math.pi:
            phis = math.pi - tmphis
            out = bTools.bktrtsafe(bTools.fnphi1, phis)
            return math.pi - out

        if tmphis >= math.pi and tmphis < math.pi + math.pi / 2:
            phis = tmphis - math.pi
            out = bTools.bktrtsafe(bTools.fnphi1, phis)
            return math.pi + out

        if tmphis >= math.pi + math.pi / 2 and tmphis <= math.pi + math.pi:
            phis = 2 * math.pi - tmphis
            out = bTools.bktrtsafe(bTools.fnphi1, phis)
            return 2 * math.pi - out

    def fnalpha(x, phis, phi2):
        """
        Calculate the value of fnalpha for the given x, phis, and phi2.

        Args:
        x: A float representing the value of x.
        phis: A float representing the value of phis.
        phi2: A float representing the value of phi2.

        Returns:
        A float representing the calculated value of fnalpha.
        """
        return math.sqrt(
            abs(math.cos(x) - math.cos(phi2) + (x - phi2) * math.sin(phis))
        )

    def abkt_m(a, b, phis):
        """computes the alpha bucket using a, b, and phis

        Args:
            a (float): start of interval
            b (float): end of interval
            phis (float): synchronous phase angle

        Returns:
            float: value of the alpha bucket
        """
        return 16 * math.sqrt(abs(b / a)) * bTools.Alpha_bkt(phis)

    def Alpha_bkt(phis):
        """computes the alpha bucket using the phase angle

        Args:
            phis (float): synchronous phase angle

        Returns:
            float: alpha bucket using gaussian
        """
        phis = bTools.generic_phis(phis)
        phi2 = bTools.phi_2_bkt(phis)
        phi1 = bTools.phi_1_bkt(phis)
        if phi1 == float("inf"):
            return float("inf")
        return (
            0.125
            * bmath.ROOT2
            * bTools.qguass_1t(bTools.fnalpha, phi1, phi2, phis, phi2)
        )

    def fnphi1bun(x, phis, phi2):
        """computes the function and derivative at x, phis, and phi2

        Args:
            x (float): value to integrate over
            phis (float): synchronous phase angle
            phi2 (float): phase angle

        Returns:
            float, float: the value at x and phis and the derivative at x and phis
        """
        fn = math.cos(x) - math.cos(phi2) + (x - phi2) * math.sin(phis)
        df = math.sin(phis) - math.sin(x)
        return fn, df

    def phi_1_bun(phis, phi2):
        """Determines where the phase angle is and computes the phase angle bunch value

        Args:
            phis (float): synchronous phase angle
            phi2 (float): phase angle

        Returns:
            float: phase angle bunch value
        """
        if phis >= 0 and phi2 < bmath.pi_2:
            if phi2 < phis or phi2 > (bmath.pi - phis):
                # return float("inf")
                return 1
            return bTools.rtsafe_l_t(
                bTools.fnphi1bun, (bmath.pi - phis), phis, phis, phi2
            )

        if phis >= bmath.pi_2 and phis <= bmath.pi:
            phis = bmath.pi - phis
            phi2 = bmath.pi - phi2
            if phi2 < phis or phi2 > (bmath.pi - phis):
                return 1
                # return float("inf")
            return bmath.pi - bTools.rtsafe_l_t(
                bTools.fnphi1bun, (-bmath.pi - phis), phis, phis, phi2
            )

        if phis >= bmath.pi and phis < (bmath.pi + bmath.pi_2):
            phis = phis - bmath.pi
            phi2 = phi2 - bmath.pi
            if phi2 < phis or phi2 > (bmath.pi - phis):
                # return float("inf")
                return 1
            return bmath.pi + bTools.rtsafe_l_t(
                bTools.fnphi1bun, (-bmath.pi - phis), phis, phis, phi2
            )

        if phis >= bmath.pi + bmath.pi_2 and phis <= (bmath.pi + bmath.pi):
            phis = 2 * bmath.pi - phis
            phi2 = 2 * bmath.pi - phi2
            if phi2 < phis or phi2 > (bmath.pi - phis):
                # return float("inf")
                return 1
            return 2 * bmath.pi - bTools.rtsafe_l_t(
                bTools.fnphi1bun, (-bmath.pi - phis), phis, phis, phi2
            )

    def alpha_bun(phis, phi2):
        """computes the generic phi2, phis and phi1 bunch to find the alpha bunch

        Args:
            phis (float): synchronous phase angle
            phi2 (float): phase angle

        Returns:
            float : alpha bunch
        """
        phi2 = bTools.generic_phi2(phis, phi2)
        phis = bTools.generic_phis(phis)
        phi1 = bTools.phi_1_bun(phis, phi2)
        return (
            0.125
            * bmath.ROOT2
            * bTools.qguass_1t(bTools.fnalpha, phi1, phi2, phis, phi2)
        )

    def phi_1_phi12(phi12, phis):
        """Computes the temporary phase angle depending on the current phase angles

        Args:
            phi12 (float): phase angle
            phis (float): synchronous phase angle

        Returns:
            float: returns the temporary value of the new phase angle
        """
        if phis < bmath.epsilon:
            return -0.5 * phi12
        if phi12 < bmath.epsilon:
            return phis
        if math.sin(0.5 * phi12) < 0.5 * phi12 * math.sin(phis):
            return float("inf")
        tmp = -0.5 * phi12 + math.asin(
            0.5 * phi12 * math.sin(phis) / math.sin(0.5 * phi12)
        )
        if (tmp + phi12) > (bmath.pi - phis):
            return float("inf")
        return tmp

    def alpha_phi12(phi12, phis):
        """computes the alpha phase angle

        Args:
            phi12 (float): phase angle
            phis (float): synchronous phase angle

        Returns:
            float: the gaussian of the phase angles
        """
        phis = bTools.generic_phis(phis)
        phi1 = bTools.phi_1_phi12(phi12, phis)
        phi2 = phi1 + phi12
        return (
            0.125
            * bmath.ROOT2
            * bTools.qguass_1t(bTools.fnalpha, phi1, phi2, phis, phi2)
        )

    def fbun(x, phi1, phi2, phis):
        """bunch function

        Args:
            x (float): current x value
            phi1 (float): phase angle 1
            phi2 (float): phase angle 2
            phis (float): synchronous phase angle
        Returns:
            float: bunch result
        """
        y1 = x - phi1
        y2 = phi2 - x
        y3 = math.cos(x) - math.cos(phi2) + (x - phi2) * math.sin(phis)
        return math.sqrt(abs(y2 * y1 / y3))

    def fbuncd(phi2, alpha, phis, fn, df):
        """Derivitive of the bunch function

        Args:
            phi2 (float): phase angle 2
            alpha (float): alpha value
            phis (float): synchronous phase angle
            fn (function): function
            df (float): derivative of the function

        Returns:
            float: derivative of bunch function using chebyshev
        """
        fn = bTools.alpha_bun(phis, phi2) - alpha
        phi1 = bTools.phi_1_bun(phis, phi2)
        df = (
            0.0625
            * bmath.ROOT2
            * (math.sin(phi2) - math.sin(phis))
            * bTools.qchebyshev_t(bTools.fbun, phi1, phi2, phis)
        )
        return fn, df

    def i_alpha_bun(alpha, phis):
        """

        Args:
            alpha (float): alpha value
            phis (float): synchronous phase angle
        Returns:
            float : synchronous phase angle or secant method approximation
        """
        phis = bTools.generic_phis(phis)
        return (
            phis
            if alpha < 1.0e-10
            else bTools.rtsafe_l_t(bTools.fbuncd, phis, (bmath.pi - phis), alpha, phis)
        )

    def fnTbun(x, phi1, phi2, phis):
        """

        Args:
            x (float): current x value on bunch/bucket plot
            phi1 (float): first phase angle
            phi2 (float): second phase angle
            phis (float): synchronous phase angle

        Returns:
            float : returns the calculated value of the bunch
        """
        y1 = x - phi1
        y2 = phi2 - x
        # print(x)
        # print(phi1)
        # print(phi2)
        # print(phis)
        if phis == float("inf"):
            phis = 1
        y3 = math.cos(x) - math.cos(phi2) + (x - phi2) * math.sin(phis)
        y3 = 1
        return math.sqrt(abs(y2 * y1 / y3))

    def T_bun(phis, phi2):
        """computes the bunch using the chebyshev approximation

        Args:
            phis (float): synchronous phase angle
            phi2 (float): second phase angle

        Returns:
            float : chebyshev approximation of the bunch
        """
        phi2 = bTools.generic_phi2(phis, phi2)
        phis = bTools.generic_phis(phis)

        p2bkt = bTools.phi_2_bkt(phis)
        p1bkt = bTools.phi_1_bkt(phis)

        # if (phi2>p2bkt || phi2<p1bkt) return HUGE_VAL; /* if outside bkt */
        if phi2 > p2bkt or phi2 < p1bkt:
            return 10e100
            # if outside bkt
        p2 = phi2 if phi2 > phis else (2 * phis - phi2)  # /* make sure p2>ps */
        phi1 = bTools.phi_1_bun(phis, p2)
        # print("phi1: "+str(phi1))
        p1 = phi1
        return bTools.qchebyshev_t(bTools.fnTbun, p1, p2, phis)

    def fnT2rfbun(x, vrf, vn, n, theta, phis, phi1, phi2):
        """computes the second rf bunch values

        Args:
            x (float): current x value on the plaot
            vrf (float): rf voltage
            vn (float): rf harmonic number
            n (int): number of bunches
            theta (float): angle value
            phis (float): synchronous phase angle
            phi1 (float): first phase angle
            phi2 (float): secong phase angle

        Returns:
            float: calculated second rf bunch value
        """
        y1 = x - phi1
        y2 = phi2 - x
        y3 = (
            vrf * (math.cos(x) - math.cos(phi2) + (x - phi2) * math.sin(phis))
            + vn / n * (math.cos(n * (x + theta)) - math.cos(n * (phi2 + theta)))
            + vn * (x - phi2) * math.sin(n * (phis + theta))
        )
        return math.sqrt(abs(y2 * y1 / y3))

    def T_2rfbun(vrf, vn, n, theta, phis, phi1, phi2):
        """Chebyshev approxmimation of the second rf bunch

        Args:
            vrf (float): rf voltage
            vn (float): rf harmonic number
            n (int): number of bunches
            theta (float): angle value
            phis (float): synchronous phase angle
            phi1 (float): first phase angle
            phi2 (float): secong phase angle

        Returns:
            float: second rf bunch value
        """
        # /* make sure phi2>phi1 */
        return bTools.qchebyshev_7t(
            bTools.fnT2rfbun, vrf, vn, n, theta, phis, phi1, phi2
        )

    def phi_2_dW(phis, A, B, dW):
        """computes the second phase angle using the wave derivative

        Args:
            phis (float): synchronous phase angle
            A (float): value for A
            B (float): value for B
            dW (float): value for dW

        Returns:
            float: value for phase angle
        """
        i = 1000
        N = 1000
        dphi = bmath.pi / N

        if A < 0:
            phi2 = bmath.pi - phis
            phi1 = phis
        else:
            phi1 = bmath.pi - phis
            phi2 = phis

        dH = 0.5 * A * dW * dW
        phi = phi2
        # /*printf("phis=%lf A=%lf B=%lf dW=%lf phi1=%lf dphi=%lf\n",phis,A,B,dW,phi1,dphi);*/
        for i in range(0, N):
            phi -= dphi
            dh = (
                B
                * (math.cos(phi) - math.cos(phis) + (phi - phis) * math.sin(phis))
                / dH
                - 1.0
            )
            if abs(dh) < 0.01:
                return phi

    ## /lib/archive/a/BktBun.c
    def U2rf(x, A, v1, vn, n, theta, phis):
        """Computes the value to use for the 2rfu graph

        Args:
            x (float): x value from plot
            A (float): value of A
            v1 (float): value of vrf
            vn (float): value of vn
            n (int): number of bunches
            theta (float): the angle theta value
            phis (float): synchronous phase angle

        Returns:
            float: 2rfu value
        """
        U = np.sign(A) * (
            v1 * (math.cos(x) - math.cos(phis) + (x - phis) * math.sin(phis))
            + vn / n * (math.cos(n * (x + theta)) - math.cos(n * (phis + theta)))
            + vn * (x - phis) * math.sin(n * (phis + theta))
        )
        return U

    def Phi_s(etas, C, rho, bdot, vrfk):
        # print(etas, C, rho, bdot, vrfk)
        tmp = C * rho * bdot / vrfk
        # print(tmp)
        if abs(tmp) > 1:
            print("Error: |C*rho*Bdot/Vrf| > 1")
            return 0

        if abs(etas) < bmath.epsilon:
            print(
                "Error: slip factor |etas| < 10^-6, in the transition regime, stay away from this region"
            )

        if vrfk > 0.0:
            out = math.asin(tmp) if etas < 0 else math.pi - math.asin(tmp)
            return out
            # Tcl_PrintDouble(interp, out, res)
            # Tcl_AppendResult(interp, res, None)
        else:
            out = 2 * math.pi - math.asin(tmp) if etas < 0 else math.pi + math.asin(tmp)
            return out
            # Tcl_PrintDouble(interp, out, res)
            # Tcl_AppendResult(interp, res, None)
            # return TCL_OK

    def Phi_2_bkt(phis, etas, C, rho, bdot, vrf):
        if phis >= 0 and phis < bmath.pi_2:
            out = bmath.pi - phis
            return out
            # Tcl_PrintDouble(interp,out,res)
            # Tcl_AppendResult(interp, res, (char*) NULL)

        if phis >= bmath.pi_2 and phis <= bmath.pi:
            out = bmath.pi - phis
            return out
            # Tcl_PrintDouble(interp,out,res)
            # Tcl_AppendResult(interp, res, (char*) NULL)
            # return TCL_OK

        if phis >= bmath.pi and phis < (bmath.pi + bmath.pi_2):
            out = 3 * bmath.pi - phis
            return out
            # Tcl_PrintDouble(interp,out,res)
            # Tcl_AppendResult(interp, res, (char*) NULL)
            # return TCL_OK

        if phis >= bmath.pi + bmath.pi_2 and phis <= (bmath.pi + bmath.pi):
            out = 3 * bmath.pi - phis
            return out
            # Tcl_PrintDouble(interp,out,res)
            # Tcl_AppendResult(interp, res, (char*) NULL)
            # return TCL_OK

    def Phi_1_bkt(phis, etas, C, rho, bdot, vrf):
        out = bTools.bktrtsafe(bTools.fnphi1, phis)

        if abs(out - float("inf")) < bmath.epsilon:
            print("Error: MAXIT is too small in Phi_1_bkt.c")

        if phis >= 0 and phis < bmath.pi_2:
            out = bTools.bktrtsafe(bTools.fnphi1, phis)
        #    Tcl_PrintDouble(interp,out,res)
        #    Tcl_AppendResult(interp, res, (char*) NULL)
        #    return TCL_OK

        if phis >= bmath.pi_2 and phis <= bmath.pi:
            phis = bmath.pi - phis
            out = bTools.bktrtsafe(bTools.fnphi1, phis)
            out = bmath.pi - out
        #    Tcl_PrintDouble(interp,out,res)
        #    Tcl_AppendResult(interp, res, (char*) NULL)
        #    return TCL_OK

        if phis >= bmath.pi and phis < (bmath.pi + bmath.pi_2):
            phis = phis - bmath.pi
            out = bTools.bktrtsafe(bTools.fnphi1, phis)
            out = bmath.pi + out
        #    Tcl_PrintDouble(interp,out,res)
        #    Tcl_AppendResult(interp, res, (char*) NULL)
        #    return TCL_OK

        if phis >= (bmath.pi + bmath.pi_2) and phis <= (bmath.pi + bmath.pi):
            phis = 2 * bmath.pi - phis
            out = bTools.bktrtsafe(bTools.fnphi1, phis)
            out = 2 * bmath.pi - out
        #    Tcl_PrintDouble(interp,out,res)
        #    Tcl_AppendResult(interp, res, (char*) NULL)
        #    return TCL_OK

    def alpha_bkt(phis, eta, C, rho, bdot, vrf):
        phis = bTools.generic_phis(phis)
        phi2 = bTools.phi_2_bkt(phis)
        phi1 = bTools.phi_1_bkt(phis)
        if phi1 == float("inf"):
            print("Error: wrong in phi_1_bkt.c")

        out = (
            0.125
            * bmath.ROOT2
            * bTools.qguass_1t(bTools.fnalpha, phi1, phi2, phis, phi2)
        )
        # Tcl_PrintDouble(interp,out,res)
        # Tcl_AppendResult(interp, res, (char*) NULL)
        # return TCL_OK

    def generate_bkt(phis, A, B):
        i = 200
        N = 200

        phi2 = bTools.phi_2_bkt(phis)
        phi1 = bTools.phi_1_bkt(phis)

        yo = math.sqrt(
            abs(
                2
                * B
                / A
                * (-2 * math.cos(phis) + (bmath.pi - 2 * phis) * math.sin(phis))
            )
        )
        wmax = yo
        wmin = -yo

        dymax = 1.2 * yo
        dymin = -dymax

        yclip = 0.2 * wmax

        if abs(phis) < bmath.epsilon:
            uplim = bmath.pi
            lowlim = -bmath.pi
            dx = (uplim - lowlim) / (N - 1)

            # fp=fopen(".gbkt.dat","w")
            # p = dvector(1,2*N)
            # w = dvector(1,2*N)
            p = []
            w = []
            for i in range(0, N):
                p.append(uplim - dx * (i - 1))
                w[i] = math.sqrt(
                    abs(
                        -2
                        * B
                        / A
                        * (
                            math.cos(p[i])
                            - math.cos(phi2)
                            + (p[i] - phi2) * math.sin(phis)
                        )
                    )
                )
                # fprintf(fp,"%lf %lf\n",p[i]*radeg,w[i])
                # /*	printf("%lf %lf\n",p[i]*radeg,w[i]);*/

            for i in range(N + 1, N * 2):
                p[i] = p[2 * N + 1 - i]
                w[i] = -w[2 * N + 1 - i]
                # fprintf(fp,"%lf %lf\n",p[i]*bmath.radeg,w[i])

            # fclose (fp)
            # free_dvector(p,1,2*N)
            # free_dvector(w,1,2*N)
            # return TCL_OK

        if abs(phis - bmath.pi) < bmath.epsilon:
            uplim = 2 * bmath.pi
            lolim = 0
            dx = (uplim - lolim) / (N - 1)

            # fp=fopen(".gbkt.dat","w")
            # p = dvector(1,2*N)
            # w = dvector(1,2*N)
            p = []
            w = []
            for i in range(0, N):
                p[i] = uplim - dx * (i - 1)
                w[i] = math.sqrt(
                    abs(
                        -2
                        * B
                        / A
                        * (
                            math.cos(p[i])
                            - math.cos(phi2)
                            + (p[i] - phi2) * math.sin(phis)
                        )
                    )
                )
                # fprintf(fp,"%lf %lf\n",p[i]*bmath.radeg,w[i])
                # /*	printf("%lf %lf\n",p[i]*radeg,w[i]);*/

            for i in range(N + 1, N * 2):
                p[i] = p[2 * N + 1 - i]
                w[i] = -w[2 * N + 1 - i]
                # fprintf(fp,"%lf %lf\n",p[i]*bmath.radeg,w[i])

        if phi2 > phi1:
            uplim = 2 * bmath.pi
            lolim = phi1
            dx = (uplim - lolim) / (N - 1)
            for i in range(0, N):
                x = uplim - dx * i
                y = math.sqrt(
                    abs(
                        -2
                        * B
                        / A
                        * (math.cos(x) - math.cos(phi2) + (x - phi2) * math.sin(phis))
                    )
                )
                if y < yclip:
                    uplim = x
                    break
            dx = (uplim - lolim) / (N - 1)
            # fp=fopen(".gbkt.dat","w")
            # p = dvector(1,2*N)
            # w = dvector(1,2*N)
            p = []
            w = []
            for i in range(0, N):
                p.append(uplim - dx * (i - 1))
                w.append(
                    math.sqrt(
                        abs(
                            -2
                            * B
                            / A
                            * (
                                math.cos(p[i])
                                - math.cos(phi2)
                                + (p[i] - phi2) * math.sin(phis)
                            )
                        )
                    )
                )
                # fprintf(fp,"%lf %lf\n",p[i]*bmath.radeg,w[i])
                # /*	printf("%lf %lf\n",p[i]*bmath.radeg,w[i]);*/

            for i in range(N + 1, N * 2):
                p[i] = p[2 * N + 1 - i]
                w[i] = -w[2 * N + 1 - i]
                # fprintf(fp,"%lf %lf\n",p[i]*bmath.radeg,w[i])
        else:
            lolim = -2 * bmath.pi
            uplim = phi2

            dx = (uplim - lolim) / (N - 1)
            for i in range(0, N):
                x = uplim - dx * (i - 1)
                y = math.sqrt(
                    abs(
                        -2
                        * B
                        / A
                        * (math.cos(x) - math.cos(phi2) + (x - phi2) * math.sin(phis))
                    )
                )
                if y < yclip:
                    lolim = x
                    uplim = phi1
                    break

            # /*	printf("%lf %lf\n",lolim,uplim);*/
            dx = (uplim - lolim) / (N - 1)

            # fp=fopen(".gbkt.dat","w")
            # p = dvector(1,2*N)
            # w = dvector(1,2*N)
            p = []
            w = []

            for i in range(0, N):
                p.append(lolim + dx * (i - 1))
                w.append(
                    (
                        abs(
                            -2
                            * B
                            / A
                            * (
                                math.cos(p[i])
                                - math.cos(phi2)
                                + (p[i] - phi2) * math.sin(phis)
                            )
                        )
                    )
                )
                # fprintf(fp,"%lf %lf\n",p[i]*radeg,w[i])
                # /*	printf("%lf %lf\n",p[i]*radeg,w[i]);*/

            for i in range(N + 1, N * 2):
                p[i] = p[2 * N + 1 - i]
                w[i] = -w[2 * N + 1 - i]
                # fprintf(fp,"%lf %lf\n",p[i]*radeg,w[i])

    def Generate_U_bkt(phis):
        i = 361
        N = 361
        dx = 4 * bmath.pi / (N - 1)
        lolim = -2 * bmath.pi
        uplim = 2 * bmath.pi
        # fp = fopen(".gUbkt.dat","w")
        eta = -1 if phis < bmath.pi_2 else 1
        for i in range(0, N):
            x = lolim + dx * i
            y = math.cos(x) - math.cos(phis) + (x - phis) * math.sin(phis)
            # fprintf(fp,"%lf %lf\n",x*radeg,y*eta);
        # close file

    def Generate_Vrf_bkt(Vrf, phis):
        i = 361
        N = 361
        amp = Vrf / Vrf
        dx = 4 * bmath.pi / (N - 1)
        lolim = -2 * bmath.pi
        uplim = 2 * bmath.pi
        # fp = fopen (".gVrfbkt.dat","w");
        eta = -1 if phis < bmath.pi_2 else 1
        for i in range(0, N):
            x = lolim + dx * i
            y = amp * (bmath.sin(x))
            # fprintf(fp,"%lf %lf\n",x*radeg,y);
        # fclose(fp);

    def Alpha_bun_phi12(phi12, phis):
        # double phi1,phi2,out,phis,phi12,phi12bkt;
        phi2 = bTools.phi_2_bkt(phis)
        phi1 = bTools.phi_1_bkt(phis)
        phi12bkt = abs(phi2 - phi1)
        phi12 = abs(phi12)
        if phi12bkt < phi12:
            print("Error: Bunch length is greater than bucket length")
        out = bTools.alpha_phi12(phi12, phis)
        return out

    def i_Alpha_bun(abun, phis):
        # double phi1,phi2,out,phis,phi12,phi12bkt,abun,abkt;
        abkt = bTools.Alpha_bkt(phis)
        if abun > abkt:
            print("Error: Bkt is smaller than bunch")

        phis = bTools.generic_phis(phis)
        # this is also not defined
        out = bTools.i_alpha_bun(abun, phis)
        # Tcl_PrintDouble(interp,out,res);
        # Tcl_AppendResult(interp, res, (char*) NULL);
        # return TCL_OK;
        return out

    def Phi_1_bun(phis, phi2):
        phis = bTools.generic_phis(phis)
        phi2 = bTools.generic_phi2(phis, phi2)
        out = bTools.phi_1_bun(phis, phi2)
        return out

    def Phi_1_phi12(phi12, phis):
        out = bTools.phi_1_phi12(phi12, phis)
        return out

    def Proper_phi(phis, phi):
        out = bTools.proper_phi(phis, phi)
        return out

    def Generic_phi2(phis, phi):
        out = bTools.generic_phi2(phis, phi)
        return out

    def Generate_bun(phis, phi1, phi2, A, B):
        i = 200
        N = 200
        uplim = phi2
        lolim = phi1
        dx = (uplim - lolim) / (N - 1)
        p = []
        w = []
        for i in range(0, 200):
            p.append(lolim + dx * (i - 1))
            w.append(
                math.sqrt(
                    abs(
                        -2
                        * B
                        / A
                        * (
                            math.cos(p[i])
                            - math.cos(phi1)
                            + (p[i] - phi1) * math.sin(phis)
                        )
                    )
                )
            )
            # fprintf(fp,"%lf %lf\n",p[i]*radeg,w[i]);

    def Period_bun(phis, phi2, A, B):
        out = math.sqrt(abs(2.0 / A / B)) * bTools.T_bun(phis, phi2)
        return out

    def RESET(a):
        N = 20
        dx = bmath.pi / (N - 1)
        p = np.zeros(N)
        for i in range(0, N):
            p[2 * i] = dx * i
            p[2 * i + 1] = a * math.sin(p[2 * i])
            p[2 * i] = p[2 * i] * 180 / bmath.pi
        # if (Blt_GraphElement(interp, ".gbkt", "BUN", 2*N, p) != TCL_OK) {
        #    Tcl_AppendResult(interp, "wrong # args: \"", argv[0]," ps", (char*) NULL);
        #    return TCL_ERROR;}

    def DrawBun(phis, phi1, phi2, A, B):
        N = 100
        uplim = phi2
        lolim = phi1
        dx = (uplim - lolim) / (N - 1)
        p = np.zeros(N)
        for i in range(0, N):
            p[2 * i] = lolim + dx * i
            p[2 * i + 1] = math.sqrt(
                abs(
                    -2
                    * B
                    / A
                    * (
                        math.cos(p[2 * i])
                        - math.cos(phi1)
                        + (p[2 * i] - phi1) * math.sin(phis)
                    )
                )
            )
            p[2 * i] *= bmath.radeg

        for i in range(N + 1, N * 2):
            p[2 * i] = p[4 * N - 2 - 2 * i]
            p[2 * i + 1] = -p[4 * N - (2 * i + 1)]

        # if (Blt_GraphElement(interp, ".gbkt", "BUN", 4*N, p) != TCL_OK) {
        # Tcl_AppendResult(interp, "wrong # args: \"", argv[0]," ps",

    #       (char*) NULL);
    # return TCL_ERROR;}

    def DrawBkt(phis, A, B):
        N = 200
        phi2 = bTools.phi_2_bkt(phis)
        phi1 = bTools.phi_1_bkt(phis)

        yo = math.sqrt(
            abs(
                2
                * B
                / A
                * (-2 * math.cos(phis) + (bmath.pi - 2 * phis) * math.sin(phis))
            )
        )
        wmax = yo
        wmin = -yo

        dymax = 1.2 * yo
        dymin = -dymax

        yclip = 0.2 * wmax

        if abs(phis) < bmath.epsilon:
            uplim = bmath.pi
            lolim = -bmath.pi
            dx = (uplim - lolim) / (N - 1)

            p = np.zeros(4 * N)
            for i in range(0, N):
                p[2 * i] = uplim - dx * i
                p[2 * i + 1] = math.sqrt(
                    abs(
                        -2
                        * B
                        / A
                        * (
                            math.cos(p[2 * i])
                            - math.cos(phi2)
                            + (p[2 * i] - phi2) * math.sin(phis)
                        )
                    )
                )
                p[2 * i] *= bmath.radeg

            for i in range(N + 1, N * 2):
                p[2 * i] = p[4 * N - 2 - 2 * i]
                p[2 * i + 1] = -p[4 * N - (2 * i + 1)]

            # if (Blt_GraphElement(interp, ".gbkt", "BKT", 4*N, p) != TCL_OK) {
            # Tcl_AppendResult(interp, "wrong # args: \"", argv[0]," ps",
            #        (char*) NULL);
            # return TCL_ERROR;}
            # /*	Tcl_VarEval (interp,".gbkt element show BKT", (char *) NULL);*/
            # free_dvector(p,0,4*N);
            # return TCL_OK;

        if abs(phis - bmath.pi) < bmath.epsilon:
            uplim = 2 * bmath.pi
            lolim = 0
            dx = (uplim - lolim) / (N - 1)

            p = np.zeros(4 * N)
            for i in range(0, N):
                p[2 * i] = uplim - dx * i
                p[2 * i + 1] = math.sqrt(
                    abs(
                        -2
                        * B
                        / A
                        * (
                            math.cos(p[2 * i])
                            - math.cos(phi2)
                            + (p[2 * i] - phi2) * math.sin(phis)
                        )
                    )
                )
            for i in range(N + 1, N * 2):
                p[2 * i] = p[4 * N - 2 - 2 * i]
                p[2 * i + 1] = p[4 * N - (2 * i + 1)]

            # if (Blt_GraphElement(interp, ".gbkt", "BKT", 4*N, p) != TCL_OK) {
            # Tcl_AppendResult(interp, "wrong # args: \"", argv[0]," ps",
            #        (char*) NULL);
            # return TCL_ERROR;}
            # /*    Tcl_VarEval (interp,".gbkt element show BKT", (char *) NULL);*/
            # free_dvector(p,0,4*N);
            # return TCL_OK;

        if phi2 > phi1:
            uplim = 2 * pi
            lolim = phi1
            dx = (uplim - lolim) / (N - 1)
            for i in range(0, N):
                x = uplim - dx * i
                y = math.sqrt(
                    abs(
                        -2
                        * B
                        / A
                        * (math.cos(x) - math.cos(phi2) + (x - phi2) * math.sin(phis))
                    )
                )
                if y < yclip:
                    uplim = x
                    break
            dx = (uplim - lolim) / (N - 1)
            p = np.zeros(N * 4)
            for i in range(0, N):
                p[2 * i] = uplim - dx * i
                p[2 * i + 1] = math.sqrt(
                    abs(
                        -2
                        * B
                        / A
                        * (
                            math.cos(p[2 * i])
                            - math.cos(phi2)
                            + (p[2 * i] - phi2) * math.sin(phis)
                        )
                    )
                )
                p[2 * i] *= bmath.radeg

            for i in range(N, N * 2):
                p[2 * i] = p[4 * N - 2 - 2 * i]
                p[2 * i + 1] = -p[4 * N - (2 * i + 1)]

            # if (Blt_GraphElement(interp, ".gbkt", "BKT", 4*N, p) != TCL_OK) {
            # Tcl_AppendResult(interp, "wrong # args: \"", argv[0]," ps",
            #        (char*) NULL);
            # return TCL_ERROR;}
            # /*    Tcl_VarEval (interp,".gbkt element show BKT", (char *) NULL);*/
            # free_dvector(p,0,4*N);
            # return TCL_OK;
        else:
            lolim = -2 * pi
            uplim = phi2  # /* from phi2 down */

            dx = (uplim - lolim) / (N - 1)
            for i in range(0, N):
                x = uplim - dx * (i - 1)
                y = math.sqrt(
                    abs(
                        -2
                        * B
                        / A
                        * (math.cos(x) - math.cos(phi2) + (x - phi2) * math.sin(phis))
                    )
                )
                if y < yclip:
                    lolim = x
                    uplim = phi1
                    break
            dx = (uplim - lolim) / (N - 1)
            p = np.zeros(4 * N)
            for i in range(0, N):
                p[2 * i] = lolim + dx * i
                p[2 * i + 1] = math.sqrt(
                    abs(
                        -2
                        * B
                        / A
                        * (
                            math.cos(p[2 * i])
                            - math.cos(phi2)
                            + (p[2 * i] - phi2) * math.sin(phis)
                        )
                    )
                )
                p[2 * i] *= bmath.radeg

            for i in range(N + 1, N * 2):
                p[2 * i] = p[4 * N - 2 - 2 * i]
                p[2 * i + 1] = -p[4 * N - (2 * i + 1)]

            # if (Blt_GraphElement(interp, ".gbkt", "BKT", 4*N, p) != TCL_OK) {
            # Tcl_AppendResult(interp, "wrong # args: \"", argv[0]," ps",
            #        (char*) NULL);
            # return TCL_ERROR;}
            # /*    Tcl_VarEval (interp,".gbkt element show BKT", (char *) NULL);*/
            # free_dvector(p,0,4*N);
            # return TCL_OK;

    def Draw2rf(v1, vn, n, theta, wphase):
        N = 540
        theta = theta / bmath.radeg
        N += 1
        # there is no defintion for wphase
        dx = wphase / (N - 1)
        p = np.zeros(2 * N)
        lolim = -bmath.pi
        uplim = bmath.pi
        for i in range(0, N):
            p[2 * i] = lolim + dx * i
            p[2 * i + 1] = v1 * math.sin(p[2 * i]) + vn * math.sin(
                n * p[2 * i] + n * theta
            )
            p[2 * i] *= bmath.radeg
        return p
        # if (Blt_GraphElement(interp, graph, element, 2*N, p) != TCL_OK) {
        # Tcl_AppendResult(interp, "wrong # args: \"", argv[0],"graph element v1....",(char*) NULL);
        # return TCL_ERROR;}
        # free_dvector(p,0,2*N);

    def Draw2rfU(A, B, v1, vn, n, theta, phis, wphase):
        # B -> e/2/pi/h
        N = 540
        theta = theta / bmath.radeg
        phis = phis / bmath.radeg
        N += 1
        dx = wphase / (N - 1)
        p = np.zeros(2 * N)
        lolim = -bmath.pi
        uplim = bmath.pi
        for i in range(0, N):
            x = lolim + dx * i
            p[2 * i] = x * bmath.radeg
            p[2 * i + 1] = bTools.U2rf(x, A, v1, vn, n, theta, phis)
        return p
        # if (Blt_GraphElement(interp, graph, element, 2*N, p) != TCL_OK) {
        # Tcl_AppendResult(interp, "wrong # args: \"",

    #       argv[0],"graph element v1.... phis",
    #       (char*) NULL);
    # return TCL_ERROR;}
    # free_dvector(p,0,2*N);

    def DrawPhis(phis, height):
        N = 2
        p = np.zeros(2 * N)
        p[0] = phis
        p[1] = -height
        p[2] = phis
        p[3] = height
        return p

    def DrawHContour(
        v1, vn, n, theta, lolim, uplim, dx, *p, phis, x, h2, y, phi2, hx, wphase
    ):
        N = 540
        k = 1
        theta = theta / bmath.radeg
        phis = phis / bmath.radeg
        phi2 = phi2 / bmath.radeg

        N += 1
        dx = wphase / (N - 1)
        p = np.zeros(4 * (N + k))
        lolim = -bmath.pi
        uplim = bmath.pi

        h2 = -(
            v1 * (math.cos(phi2) - math.cos(phis) + (phi2 - phis) * math.sin(phis))
            + vn / n * (math.cos(n * (phi2 + theta)) - math.cos(n * (phis + theta)))
            + vn * (phi2 - phis) * math.sin(n * (phis + theta))
        )
        for i in range(0, N + k):
            x = lolim + dx * i
            p[2 * i] = x * bmath.radeg
            hx = -(
                v1 * (math.cos(x) - math.cos(phis) + (x - phis) * math.sin(phis))
                + vn / n * (math.cos(n * (x + theta)) - math.cos(n * (phis + theta)))
                + vn * (x - phis) * math.sin(n * (phis + theta))
            )

            y = h2 - hx
            p[2 * i + 1] = 10 * math.sqrt(abs(y)) if y > 0.0 else 0

        for i in range(0, N + k):
            p[2 * i] = p[4 * (N + k) - 2 - 2 * i]
            p[2 * i + 1] = -p[4 * (N + k) - (2 * i + 1)]

    def DrawHLines(
        v1, vn, n, theta, lolim, uplim, dx, *p, phis, x, h2, y, phi2, hx, wphase
    ):
        # double v1,vn,n,theta,lolim,uplim,dx,*p,phis,x,h2,y,phi2,hx;
        N = 540
        k = 1

        theta = theta / bmath.radeg
        phis = phis / bmath.radeg
        phi2 = phi2 / bmath.radeg

        N += 1
        dx = wphase / (N - 1)
        p = np.zeros(4 * (N + k))
        lolim = -bmath.pi
        uplim = bmath.pi

        h2 = -(
            v1 * (math.cos(phi2) - math.cos(phis) + (phi2 - phis) * math.sin(phis))
            + vn / n * (math.cos(n * (phi2 + theta)) - math.cos(n * (phis + theta)))
            + vn * (phi2 - phis) * math.sin(n * (phis + theta))
        )

        for i in range(0, N + k):
            x = lolim + dx * i
            p[2 * i] = x * bmath.radeg
            hx = -(
                v1 * (math.cos(x) - math.cos(phis) + (x - phis) * math.sin(phis))
                + vn / n * (math.cos(n * (x + theta)) - math.cos(n * (phis + theta)))
                + vn * (x - phis) * math.sin(n * (phis + theta))
            )
            y = h2 - hx
            p[2 * i + 1] = 10 * math.sqrt(abs(y)) if y > 0.0 else 0

        for i in range(0, 2 * (N + k)):
            p[2 * i] = p[4 * (N + k) - 2 - 2 * i]
            p[2 * i + 1] = -p[4 * (N + k) - (2 * i + 1)]

    def Draw2rfSep(n, wphase, theta, phis, A, v1, vn, B):
        N = 540
        k = 1
        nn = 0
        N = N * n + 1
        dx = wphase / (N - 1)
        lolim = -pi
        uplim = pi

        theta = theta / bmath.radeg
        phis = phis / bmath.radeg

        phix2 = np.zeros(int(n))
        ymax = np.zeros(int(n))
        yp = np.zeros(int(n))
        y = np.zeros(int(N))
        y[0] = bTools.U2rf(lolim - dx, A, v1, vn, n, theta, phis)
        all_Pdata = []

        for i in range(1, N):  # find the U curve
            y[i] = bTools.U2rf(lolim + dx * (i - 1), A, v1, vn, n, theta, phis)

        # print(N)
        for i in range(0, N - 1):  # find all local max
            nn = 0
            if (y[i] >= y[i - 1]) and (y[i] > y[i + 1]):
                phix2[nn] = lolim + (i - 1.0) * dx
                ymax[nn] = y[i]
                yp[nn] = i
                nn += 1

        nn -= 1  # /* nn=n-1 */
        max = ymax[0]
        for j in range(1, nn + 1):
            if ymax[j] > max:
                max = ymax[j]

        j = 0
        if (
            ymax[0] == max and nn != 0 and y[1] < max and y[N - 1] < max
        ):  # /* don't miss a well */
            j = 1

        if ymax[nn] == max and nn != 0 and y[1] < max and y[N - 1] < max:
            nn -= 1

        p = [0] * int(4 * (N + k))
        for ji in range(j, nn):  # /* see above if */
            for i in range(0, N + k):  # /* setting zero line*/
                p[2 * i] = (
                    lolim + dx * (i - 1)
                ) * bmath.radeg  # /* in degrees for blt */
                p[2 * i + 1] = 0

            ib = 0
            ie = N + k  # /* locate turning points ib,ie */
            # /*	ib= (ymax[j]>y[0]) ? yp[0] : 0;
            # ie= (ymax[j]<y[N])? yp[nn] : (N+k);

            for ibe in range(0, j):
                if ymax[ji] < ymax[ibe]:
                    ib = yp[ibe]

            for ibe in range(ji, nn, -1):
                if ymax[ji] <= ymax[ibe]:
                    ie = yp[ibe]

            phi2 = phix2[ji]
            h2 = bTools.U2rf(phi2, A, v1, vn, n, theta, phis)
            # /* p[2*(N+k)] fake*/
            for i in range(ib, ie):
                x = p[2 * i] / bmath.radeg
                # /* in radians */
                hx = bTools.U2rf(x, A, v1, vn, n, theta, phis)
                w = h2 - hx
                if w > 0:
                    p[2 * i + 1] = math.sqrt(abs(w * 2 * B / A))
            p[2 * int(yp[j]) + 1] = 0

            # /* need to clean up the curve a bit */
            # /* from left */
            for i in range(0, N + k):
                if p[2 * i + 1] == 0:
                    i = N + k
                else:
                    p[2 * i + 1] = 0

            # /* from right */
            for i in range(N + k - 1, 0, -1):
                if p[2 * i + 1] == 0:
                    i = 0
                p[2 * i + 1] = 0

            # /* areas under the top bkt? */
            # /* very tricky */

            ib = 0
            bkt = 0  # /* reuse variable ib */
            for i in range(1, N + 1):
                if p[2 * (i - 1) + 1] == 0 and p[2 * i + 1] != 0:
                    ib += 1
                    # /*printf ("non-empty bkt begin = %d \n",i);*/

                bkt += p[2 * i + 1]
                if p[2 * i + 1] != 0 and p[2 * (i + 1) + 1] == 0:
                    # /*		printf ("non-empty bkt end = %d and bkt = %lf\n",i,bkt);*/
                    bkt = 0

            # /* duplicate bottom bkt for plotting*/
            for i in range(N + k, 2 * (N + k)):
                p[2 * i] = p[4 * (N + k) - 2 - 2 * i]
                p[2 * i + 1] = -p[4 * (N + k) - (2 * i + 1)]

            # all_Pdata.append(p)
        # print(all_Pdata)
        # print(phix2)
        # print(ymax)
        # print(yp)
        # print(y)
        # return all_Pdata
        return p, y

    def Find2rfSep(n, wphase, theta, phis, A, v1, vn):
        # double v1,vn,n,theta,lolim,uplim,dx,*p,phis,x,h2,w,phi2,hx,*y,*phix2;
        # double A;
        N = 540
        k = 1
        nn = 0
        N = N * n + 1
        dx = wphase / (N - 1)
        lolim = -bmath.pi
        uplim = bmath.pi

        theta = theta / bmath.radeg
        phis = phis / bmath.radeg

        dx = wphase / (N - 1)
        phix2 = np.zeros(n)
        y = np.zeros(N)
        y[0] = bTools.U2rf(lolim - dx, A, v1, vn, n, theta, phis)

        for i in range(1, N):  # /* find the U curve */
            y[i] = bTools.U2rf(lolim + dx * (i - 1), A, v1, vn, n, theta, phis)

        for i in range(1, N):  # /* find all local max */
            if (y[i] >= y[i - 1]) and (y[i] > y[i + 1]):
                phix2[nn] = lolim + (i - 1) * dx
                nn += 1
        nn -= 1  # /* nn=n-1 */

    def BKT2rf(wphase, n, theta, phis, A, v1, vn, B, nu):
        # double v1,vn,n,theta,lolim,uplim,dx,*p,phis,x,h2,w,phi2,A,B;
        # double *ymax,hx,*y,*phix2,max,bkt,Nu;
        # int i,N=540,k=1,nn,j,*yp,ib,ie,ibe;
        N = 540
        k = 1
        nn = 0
        N = N * n + 1
        dx = wphase / (N - 1)
        lolim = -bmath.pi
        uplim = bmath.pi

        theta = theta / bmath.radeg
        phis = phis / bmath.radeg

        phix2 = np.zeros(n)
        ymax = np.zeros(n)
        yp = np.zeros(n)
        y = np.zeros(N)
        y[0] = bTools.U2rf(lolim - dx, A, v1, vn, n, theta, phis)

        for i in range(1, N):  # /* find the U curve */
            y[i] = bTools.U2rf(lolim + dx * (i - 1), A, v1, vn, n, theta, phis)

        for i in range(1, N - 1):  # /* find all local max */
            if (y[i] >= y[i - 1]) and (y[i] > y[i + 1]):
                phix2[nn] = lolim + (i - 1.0) * dx
                ymax[nn] = y[i]
                yp[nn] = i
                nn += 1
        nn -= 1  # /* nn=n-1 */

        max = ymax[0]
        for j in range(0, nn):
            if ymax[j] > max:
                max = ymax[j]

        j = 0
        if (
            ymax[0] == max and nn != 0 and y[1] < max and y[N] < max
        ):  # /* don't miss a well */
            j = 1

        if ymax[nn] == max and nn != 0 and y[1] < max and y[N] < max:
            nn -= 1

        for j in range(j, nn):  # /* see above if */
            p = np.zeros(4 * (N + k))
            for i in range(N + k):  # /* setting zero line*/
                p[2 * i] = (lolim + dx * (i - 1)) * bmath.radeg
                # /* in degrees for blt */
                p[2 * i + 1] = 0

            ib = 0
            ie = N + k  # /* locate turning points ib,ie */
            # /*	ib= (ymax[j]>y[0]) ? yp[0] : 0;
            # ie= (ymax[j]<y[N])? yp[nn] : (N+k);
            # */

            for ibe in range(0, j):
                if ymax[j] < ymax[ibe]:
                    ib = yp[ibe]

            for ibe in range(nn, j, -1):
                if ymax[j] < ymax[ibe]:
                    ie = yp[ibe]

            phi2 = phix2[j]
            h2 = bTools.U2rf(phi2, A, v1, vn, n, theta, phis)
            for i in range(ib, ie):  # /* p[2*(N+k)] fake*/
                x = p[2 * i] / bmath.radeg
                # /* in radians */
                hx = bTools.U2rf(x, A, v1, vn, n, theta, phis)
                w = h2 - hx
                if w > 0:
                    p[2 * i + 1] = math.sqrt(w)
            p[2 * yp[j] + 1] = 0

            # /* need to clean up the curve a bit */
            # /* from left */
            for i in range(0, N + k):
                if p[2 * i + 1] == 0:
                    i = N + k
                else:
                    p[2 * i + 1] = 0

            # /* from right */
            for i in range(N + k - 1, 0, -1):
                if p[2 * i + 1] == 0:
                    i = 0
                p[2 * i + 1] = 0
            # /* areas under the top bkt? */
            # /* very tricky */

            ib = 0
            bkt = 0
            # /* reuse variable ib */
            for i in range(1, N):
                if p[2 * (i - 1) + 1] == 0 and p[2 * i + 1] != 0:
                    ib += 1  # /*printf ("non-empty bkt begin = %d \n",i);*/

                bkt += p[2 * i + 1]
                if p[2 * i + 1] != 0 and p[2 * (i + 1) + 1] == 0:
                    bkt *= 2 * dx * math.sqrt(abs(2 * B / A)) / nu
                    bkt = 0
        return y[0]

    def Draw2rfHcontour(theta, phis, phi2, n, wphase, A, v1, vn, B):
        N = 540
        k = 1

        theta = theta / bmath.radeg
        phis = phis / bmath.radeg
        phi2 = phi2 / bmath.radeg

        N = N * n + 1
        dx = wphase / (N - 1)
        p = np.zeros(4 * (int(N) + k))
        lolim = -bmath.pi
        uplim = -bmath.pi

        for i in range(0, N + k):  # /* setting zero line*/
            p[2 * i] = (lolim + dx * (i - 1)) * bmath.radeg
            # /* in degrees for blt */
            p[2 * i + 1] = 0

        h2 = bTools.U2rf(phi2, A, v1, vn, n, theta, phis)

        for i in range(0, N + k):  # /* p[2*(N+k)] fake*/
            x = p[2 * i] / bmath.radeg
            hx = bTools.U2rf(x, A, v1, vn, n, theta, phis)
            y = h2 - hx
            if y > 0:
                p[2 * i + 1] = math.sqrt(abs(y * 2 * B / A))

        # /* massaging the curve */
        # /* need to clean up the curve a bit */
        # /* from left */
        for i in range(0, N + k):
            if p[2 * i + 1] == 0:
                i = N + k
            else:
                p[2 * i + 1] = 0

        # /* from right */
        for i in range(N + k - 1, i, -1):
            if p[2 * i + 1] == 0:
                i = 0
            p[2 * i + 1] = 0

        for i in range(N + k, 2 * (N + k)):
            p[2 * i] = p[4 * (N + k) - 2 - 2 * i]
            p[2 * i + 1] = -p[4 * (N + k) - (2 * i + 1)]
        return p

    def BUN2rf(theta, phis, phi2, n, wphase, A, vrf, vn, B, Nu):
        N = 540
        k = 1
        theta = theta / bmath.radeg
        phis = phis / bmath.radeg
        phi2 = phi2 / bmath.radeg

        N = N * n + 1
        dx = wphase / (N - 1)
        W = np.zeros(4 * (N + k))
        lolim = -bmath.pi
        uplim = bmath.pi

        for i in range(0, N + k):  # /* setting zero line*/
            W[2 * i] = (lolim + dx * (i - 1)) * bmath.radeg  # /* in degrees for blt */
            W[2 * i + 1] = 0

        h2 = bTools.U2rf(phi2, A, vrf, vn, n, theta, phis)
        for i in range(0, N + k):  # /* W[2*(N+k)] fake*/
            x = W[2 * i] / bmath.radeg
            hx = bTools.U2rf(x, A, vrf, vn, n, theta, phis)
            y = h2 - hx
            if y > 0:
                W[2 * i + 1] = math.sqrt(abs(y * 2 * B / A))

        # /* massaging the curve */
        # /* need to clean up the curve a bit */
        # /* from left */
        for i in range(0, N + k):
            if W[2 * i + 1] == 0:
                i = N + k
            else:
                W[2 * i + 1] = 0
        # /* from right */
        for i in range(N + k, 0, -1):
            if W[2 * i + 1] == 0:
                i = 0
            W[2 * i + 1] = 0

        # /* areas under the top bun? */
        # /* very tricky */
        ib = 0
        bun = 0
        Tnu = 0  # /* reuse variable ib */
        for i in range(1, N):
            if W[2 * (i - 1) + 1] == 0 and W[2 * i + 1] != 0:
                ib += 1  # /*printf ("non-empty bun begin = %d \n",i);*/
            bun += W[2 * i + 1]
            if W[2 * i + 1] != 0 and W[2 * (i + 1) + 1] == 0:
                bun *= 2 * dx / Nu
                # /* W is already normalized */
                bun = 0

    def fnu2rf(theta, phis, phi2, n, wphase, A, vrf, vn, B):
        N = 540
        k = 1
        theta = theta / bmath.radeg
        phis = phis / bmath.radeg
        phi2 = phi2 / bmath.radeg

        N = N * n + 1
        dx = wphase / (N - 1)
        W = np.zeros(4 * (N + k))
        lolim = -bmath.pi
        uplim = bmath.pi

        for i in range(0, N + k):  # /* setting zero line*/
            W[2 * i] = (lolim + dx * (i - 1)) * bmath.radeg  # /* in degrees for blt */
            W[2 * i + 1] = 0

        h2 = bTools.U2rf(phi2, A, vrf, vn, n, theta, phis)

        for i in range(0, N + k):  # /* W[2*(N+k)] fake*/
            x = W[2 * i] / bmath.radeg
            hx = bTools.U2rf(x, A, vrf, vn, n, theta, phis)
            y = h2 - hx
            if y > 0:
                W[2 * i + 1] = math.sqrt(abs(y * 2 * B / A))

        # /* massaging the curve */
        # /* need to clean up the curve a bit */
        # /* from left */
        for i in range(0, N + k):
            if W[2 * i + 1] == 0:
                i = N + k
            else:
                W[2 * i + 1] = 0

        # /* from right */
        for i in range(N + k - 1, 0, -1):
            if W[2 * i + 1] == 0:
                i = 0
            W[2 * i + 1] = 0

        # /* areas under the top bun? */
        # /* very tricky */

        ib = 0
        Tnu = 0  # /* reuse variable ib */
        for i in range(1, N):
            if W[2 * (i - 1) + 1] == 0 and W[2 * i + 1] != 0:
                ib += 1  # /*printf ("non-empty bun begin = %d \n",i);*/
                phi1 = W[2 * i] / bmath.radeg

            if W[2 * i + 1] != 0 and W[2 * (i + 1) + 1] == 0:
                # /*printf ("phi1= %lfphi2= %lf W[2*i]= %lf\n",phi1*radeg,phi2*radeg,W[2*i]);*/
                phi2 = W[2 * i] / bmath.radeg
                Tnu = math.sqrt(abs(2 / A / B)) * bTools.T_2rfbun(
                    vrf, vn, n, theta, phis, phi1, phi2
                )
                Tnu = 0

    def Phi_2_dW(phis, A, B, dW):
        phi = bTools.phi_2_dW(phis, A, B, dW)
