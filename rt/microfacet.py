import math
import vec3

#
# Microfacet
#

def Clamp(w, min, max):
    if w < min:
        return min
    if w > max:
        return max
    return w

def CosTheta(w):
    return w.z

def Cos2Theta(w):
    return w.z * w.z
    
def Sin2Theta(w):
    return max(0, 1.0-Cos2Theta(w))

def SinTheta(w):
    return math.sqrt(Sin2Theta(w))
def TanTheta(w):
    return SinTheta(w)/CosTheta(w)
def Tan2Theta(w):
    return Sin2Theta(w)/Cos2Theta(w)

def CosPhi(w):
    sinTheta = SinTheta(w)
    return 1 if sinTheta == 0 else Clamp(w.x/sinTheta, -1, 1)

def SinPhi(w):
    sinTheta = SinTheta(w)
    return 0 if sinTheta == 0 else  Clamp(w.y/sinTheta, -1, 1)

def Cos2Phi(w):
    return CosPhi(w) * CosPhi(w)

def Sin2Phi(w):
    return SinPhi(w) * SinPhi(w)

def CosDPhi(wa, wb):
    return Clamp(
        (wa.x * wb.x + wa.y * wb.y) / math.sqrt((wa.x * wa.x + wa.y * wa.y) *
        (wb.x * wb.x + wb.y * wb.y)),
        -1, 1)


def ErfInv(x): 
    x = Clamp(x, -.99999, .99999);
    w = -math.log((1.0 - x) * (1.0 + x));
    if w < 5:
        w = w - 2.5;
        p = 2.81022636e-08;
        p = 3.43273939e-07 + p * w;
        p = -3.5233877e-06 + p * w;
        p = -4.39150654e-06 + p * w;
        p = 0.00021858087 + p * w;
        p = -0.00125372503 + p * w;
        p = -0.00417768164 + p * w;
        p = 0.246640727 + p * w;
        p = 1.50140941 + p * w;
    else: 
        w = math.sqrt(w) - 3;
        p = -0.000200214257;
        p = 0.000100950558 + p * w;
        p = 0.00134934322 + p * w;
        p = -0.00367342844 + p * w;
        p = 0.00573950773 + p * w;
        p = -0.0076224613 + p * w;
        p = 0.00943887047 + p * w;
        p = 1.00167406 + p * w;
        p = 2.83297682 + p * w;
    return p * x;

def SphericalDirection (sinTheta, cosTheta, phi):
    return vec3.Vec3(sinTheta*math.cos(phi), sinTheta * math.sin(phi), cosTheta)


class Microfacet:
    def __init__(self, alpha_x, alpha_y):
        self.alpha_x = alpha_x
        self.alpha_y = alpha_y

    def BeckmannSample11(cosThetaI, U1, U2):
        if cosThetaI > .9999:
            r = math.sqrt(-math.log(1.0-U1))
            sinPhi = math.sin(2 * math.pi * U2)
            cosPhi = math.cos(2 * math.pi * U2)
            slope_x = r * cosPhi
            slope_y = r * sinPhi
            return (slope_x, slope_y)

        sinThetaI = math.sqrt(max(0.0, 1.0 - cosThetaI * cosThetaI))
        tanThetaI = sinThetaI/cosThetaI
        cotThetaI = 1.0/tanThetaI

        a = -1
        c = math.erf(cotThetaI)
        sample_x = max(U1,1e-6)
        thetaI = math.acos(cosThetaI)
        fit = 1.0 + thetaI * (-.876 + thetaI * (.4265 - .0594 * thetaI))
        b = c + (1.0+c) * math.pow(1.0-sample_x, fit)

        SQRT_PI_INV = 1.0/math.sqrt(math.pi)
        normalization = 1.0/(1.0 + c + SQRT_PI_INV * tanThetaI * math.exp(-cotThetaI * cotThetaI))

        it = 0
        for it in range(10):
            if not (b >= a and b <= c):
                b = .5 * (a + c)

            invErf = ErfInv(b)
            value = normalization * (1.0 + b + SQRT_PI_INV * tanThetaI * math.exp(-invErf * invErf)) - sample_x
            derivative = normalization * (1.0 - invErf * tanThetaI)

            if math.fabs(value) < 1e-5:
                break

            if value > 0:
                c = b
            else:
                a = b
            b -= value/derivative

        slope_x = ErfInv(b)
        slope_y = ErfInv(2.0 * max(U2, 1e-6) - 1.0)
        return (slope_x, slope_y)

    def BeckmannSampleVisible(self, wi, U1, U2): 
        wiStretch = vec3.Vec3(self.alpha_x * wi.x, self.alpha_y * wi.y, wi.z).norm()
        (slope_x, slope_y) = BeckmannSample11(CosTheta(wiStretch), U1, U2)
        tmp = CosPhi(wiStretch) * slope_x -SinPhi(wiStretch) * slope_y
        slope_y = SinPhi(wiStretch) * slope_x + CosPhi(wiStretch) * slope_y
        slope_x = tmp

        slope_x = alpha_x * slope_x
        slope_y = alpha_y * slope_y

        return vec3.Vec3(-slope_x, -slope_y, 1).norm()

    def BeckmannSample(self, wi, U1, U2):
        if self.alpha_x == self.alpha_y:
            logSample = math.log(1.0-U1)
            tan2Theta = -self.alpha_x * self.alpha_x * logSample
            phi = U2 * 2.0 * math.pi

        cosTheta = 1.0/math.sqrt(1.0 + tan2Theta)
        sinTheta = math.sqrt(max(0.0, 1.0 - cosTheta * cosTheta))
        wh = SphericalDirection(sinTheta, cosTheta, phi)
        if wi.z * wh.z < 0:
            wh = -wh
            cosTheta = -cosTheta
        return (wh, cosTheta,  phi)


    def Sample_wh(self, wo, u):
        (wh, cosTheta, phi) = self.BeckmannSample(wo, u[0], u[1])
        return (wh, cosTheta, phi)
    
    def D(self, wh):
        tan2Theta = Tan2Theta(wh)
        if math.isinf(tan2Theta):
            return 0
        cos4Theta = Cos2Theta(wh)
        cos4Theta *= cos4Theta
        return math.exp(-tan2Theta * (Cos2Phi(wh) /(self.alpha_x * self.alpha_y) + 
                        Sin2Phi(wh)/(self.alpha_y * self.alpha_y)))/ (math.pi + self.alpha_x * self.alpha_y * cos4Theta)


    def GLambda(self, w):
        absTanTheta = math.fabs(TanTheta(w))
        if math.isinf(absTanTheta):
            return 0.0
        alpha = math.sqrt(Cos2Phi(w) * self.alpha_x * self.alpha_x + Sin2Phi(w) * self.alpha_y * self.alpha_y)
        if alpha == 0 or absTanTheta == 0:
            return 0.0
        a = 1.0/(alpha * absTanTheta)
        if a > 1.6:
            return 0.0
        return (1 - 1.259 * a + 0.396 * a * a) / (3.535 * a + 2.181 * a * a)

    def G1(self, w):
        return 1.0/(1.0+self.GLambda(w))

    def G(self, wo, wi):
        return self.G1(wo) * self.G1(wi)

    def Pdf(self, wh):
        return self.D(wh) * math.fabs(CosTheta(wh))
