import microfacet
import zipinPaper
import zipin
import math
import vec3
import ray
import weight
import collections
import random
import matrix

debug = 0

def dprint(s):
    if debug == 1:
        print(s)

def reflect(wo, n):
    return -wo + n * vec3.dot(wo, n) * 2.0

def checkCorrectO(O, I, H):
    midOI = (O + I).norm()
    diff = H - midOI
    dis = diff.dot(diff)
    if  dis > .00001:
        print (O, I, midOI, H)
    return dis < .00001


##################################################################
#       compute theta (angle between w and N)                   ##
##################################################################

def computeTheta(w, phi = None):
    if phi == None:
        phi = math.atan2(w.y, w.x)
    cosTheta = w.z
    cosPhi = math.cos(phi)
    sinPhi = math.sin(phi)
    if math.fabs(cosPhi) > microfacet.SMALLFLOAT:
        sinTheta = w.x/cosPhi
    else:
        if math.fabs(sinPhi) > microfacet.SMALLFLOAT:
            sinTheta = w.x/sinPhi
        else:
            #raise exception
            return 0
    theta = math.atan2(sinTheta, cosTheta)
    return theta

def makeWi(theta_i, groovePhi, wirY):
    thetaI = (math.radians(theta_i));
    if thetaI < 0:
        thetaI = math.fabs(thetaI) *  math.copysign(1, groovePhi)
    else: 
        thetaI = -math.fabs(thetaI) * math.copysign(1, groovePhi)
    sinThetaI = math.sin(thetaI)
    cosThetaI = math.cos(thetaI)
    wip = microfacet.SphericalDirection(sinThetaI, cosThetaI, 0)
    #thetaIP = computeTheta(wip, 0)
    #assert(math.fabs(thetaIP-thetaI) < .00001)
    wir = wip
    len = 1-wirY * wirY
    wir *= math.sqrt(len)
    wir.y = wirY
    thetaIP = computeTheta(wir)
    return (wip, wir, thetaIP)
     


#TODO: reorganize so Brdf is base class and beckmann Brdf becomes a subclass of Brdf
#
class Brdf:
    def __init__(self, alpha_x, alpha_y):
        self.microfacet = microfacet.Beckmann(alpha_x, alpha_y)
   
    def Pdf(self, wo, wh):
        normFactor = 4.0 * vec3.dot(wo, wh)
        if normFactor < microfacet.SMALLFLOAT:
            return 0
        pdf = self.microfacet.Pdf(wh)
        return pdf/normFactor


    #given wo, return sampled wi

    def MicrofacetValue(self, wo, wi, wh, GMode = 1):
        cosThetaI = microfacet.CosTheta(wi)
        cosThetaO = microfacet.CosTheta(wo)
        if cosThetaI == 0 or cosThetaO == 0:
            return 0
        if wh.x==0 and wh.y ==0 and wh.z == 0:
            return 0
        value = self.microfacet.D(wh)
        value /= 4.0 * cosThetaO
        return  value * self.microfacet.G(wo, wi, wh, GMode)
        
    def Sample(self, wo, u, GMode = 1):
        wo = wo.norm()
        if wo.z == 0:
            return (0, 0, None)
        (wh, cosTheta, phi) = self.microfacet.Sample_wh(wo, u)
        if wh.x == 0 and wh.y == 0 and wh.z == 0:
            return (0, 0, None)
        wi = reflect(wo, wh)
        theta = computeTheta(wi)
        (value, pdf) = self.Eval(wo, wi, GMode)
        return (value, pdf, wi, theta)
    

    def Eval(self, wo, wi, GMode = 1):
        wo = wo.norm()
        wh = (wo + wi).norm()
        value = self.MicrofacetValue(wo, wi, wh, GMode)
        pdf = self.Pdf(wo, wh)
        return (value, pdf)



ANGLEERROR = 1e-6 
MAXBOUNCE = 10 
MICROFACETTHRESH = .000001
MAXITER = 1


##################################################################
#       compute phi:  we will not use the result of this
#       function for now as we are assuming phi == 0 for our    ##
#       current tests.                                          ##
##################################################################
def computePhi(w):
    cosTheta = w.z
    sinTheta = math.sqrt(1.0 - w.z * w.z)
    if math.fabs(sinTheta) < microfacet.SMALLFLOAT:
        return (0, 0)
    sinPhi = w.y/sinTheta
    cosPhi = w.x/sinTheta
    phi1 = math.atan2(sinPhi, cosPhi)
    phi2 = math.atan2(-sinPhi, -cosPhi)
    return (phi1, phi2)
   
def computeZipinNormal(grooveTheta, side, wo):
    #if grooveTheta > math.pi * .5:
    #    print('error!')
    n = vec3.Vec3(math.cos(grooveTheta), 0, math.sin(grooveTheta))
    #feng comment please why the sign works the way it does
    #will this have impact when wo.y != 0 
    n.x *= math.copysign(1, wo.x) 
    if side == 'right':
        n.x *= -1
    return n
 

class ZipinBrdf(Brdf):
    def __init__(self, alpha_x, alpha_y, zipinVersion):
        self.microfacet = microfacet.Beckmann(alpha_x, alpha_y)
        print('zipinVersion: ', zipinVersion)
        self.version = zipinVersion
        if zipinVersion == 'Feng':
            self.zipin = zipinPaper.zipin
        else:
            self.zipin = zipin.zipin

     
    def zipinRatio(self, wo, wh, prob):
        denom = math.fabs(vec3.dot(wo, wh))
        if denom < microfacet.SMALLFLOAT:
            return 1
        nom = math.fabs(wo.z * wh.z * 2.0)
        if nom < microfacet.SMALLFLOAT:
            return 0 
        gr = nom/denom * prob
        if gr > 1:
            gr = 1
        return gr 
        

    #TODO:  figure out the correct normalization factors
    #       and angle orientations for grooveTheta
    #       groovePhi
    #       and final grooveChi
    def Sample(self, wo, u, maxBounce = MAXBOUNCE):
        wo = wo.norm()
        if wo.z <=  microfacet.SMALLFLOAT:
            return (0, 0, None)
        
        for i in range(MAXITER): 
            (wh, cosThetaH, phiH) = self.microfacet.Sample_wh(wo, u)
            if wh.z < microfacet.SMALLFLOAT:
                return (0, 0, None, None)
        
            wi = reflect(wo, wh)
            (value, pdf) = self.Eval(wo, wi, maxBounce, 1)
            thetaI  = computeTheta(wi)
            return (value, pdf, wi, thetaI)

            '''
          
            #print('phiH: ' + repr(math.degrees(phiH))) 
            rotation = matrix.Rotate(wh)
            whr = rotation.rotate(wh)
            wop = wor = rotation.rotate(wo)

            #print (whr, wh)
            #whr = wh
            #wop = wor = wo
            wop.y = 0  
            wop = wop.norm()

            groovePhi = computeTheta(wop, 0)
            grooveAlpha = computeTheta(whr, 0)
            grooveTheta = .5 * math.pi - math.fabs(grooveAlpha)
            
            sampleGroovePhi = math.degrees(math.fabs(groovePhi))
            sampleGrooveTheta = math.degrees(grooveTheta)  
            hits = self.zipin(sampleGrooveTheta, sampleGroovePhi)
            (sample, prob) = weight.weightedRandomChoice(hits, maxBounce)
            if sample != None:
                (wip, wir, thetaI) = makeWi(theta_i=sample[0], groovePhi=groovePhi, wirY = -wor.y)
                (value, pdf) = self.Eval(wor, wir, maxBounce, 1)
                wi = rotation.inverse_rotate(wir)
                #wi = wir
                return (value, pdf, wi, thetaI)
            '''
        return (0, 0, None, None)



    #################################################################
    ######  Eval Helpers                            #################
    #       needs work, don't test yet                              #
    #################################################################
    def EvalNearHits(self, wo, wi, phi, chi, n, wop, wip):
        xchi = chi 
        if (n+1)%2 == 1:
            xchi *= -1
        theta = (math.pi + phi - xchi) *.5/n
        #print('right ' + repr(n))
        #print (math.degrees(theta))
        if theta <= phi or theta >= .5 * math.pi:
            return ((False, 0, 0), theta)
        return (self.EvalHits(wo, wi, phi, chi, theta, n, 'right', wop, wip), theta)

    def EvalFarHits(self, wo, wi, phi, chi, n, wop, wip):
        xchi = chi
        if n%2 == 1:
            xchi *= -1
        theta = (math.pi -phi -xchi) * .5/n
        #print('left ' + repr(n))
        #print (math.degrees(theta))
        if theta < microfacet.SMALLFLOAT or theta >=.5 * math.pi:
            return ((False, 0, 0), theta)
        return (self.EvalHits(wo, wi, phi, chi, theta, n, 'left', wop, wip), theta)

    def EvalHits(self, wo, wi, groovePhi, grooveChi,  grooveTheta, n, side, wop, wip):
        
        #0< grooveTheta < pi * .5       
        wh = computeZipinNormal(grooveTheta, side, wop)
        microfacetValue = self.MicrofacetValue(wo, wi, wh, GMode = 0)
        microfacetPdf = self.Pdf(wo, wh)
        hits = self.zipin(math.degrees(grooveTheta), math.degrees(groovePhi))
        #print(math.degrees(grooveTheta), math.degrees(groovePhi), math.degrees(grooveChi))
        #print(wh, n, side)
        #print(microfacetValue, microfacetPdf)
        #print(hits)

        if microfacetValue < 1e-6:
            return (False, 0, 0)
       
        value = 0;
        pdf = 0;
        for hit in hits:
            if hit[3] != side:
                continue
            hitChi = math.radians(hit[0])
            diff = math.fabs(hitChi - grooveChi)
            if diff < ANGLEERROR and hit[2] == n:
                hitScale =  hit[4] 
                if n == 1: 
                    blin = self.microfacet.GBlin(wo, wi, wh) 
                    if math.fabs(blin - hit[4])> .00001 or checkCorrectO(wo, wi, wh)==0:
                        print('wo, wi, wh')
                        print(wo, wi, wh)
                        print('gTheta, gPhi, gChi')
                        print(math.degrees(grooveTheta), math.degrees(groovePhi), math.degrees(grooveChi))
                        print('hitChi, gChi')
                        print(hitChi, grooveChi)
                        print('brdf, pdf of wh')
                        print(microfacetValue, microfacetPdf)
                        print('bounce count ' + repr(n))
                        GR = self.zipinRatio(wop, wh, hit[1])
                        print('g terms')
                        print(blin, hit[4], GR)
                        print('hits')
                        print(hits)
                value += microfacetValue *hitScale
                pdf += microfacetPdf * hitScale
                return (True, value, pdf)
        return (False, 0, 0)  
        

    def Eval(self, wo, wi, maxBounce = MAXBOUNCE, minBounce = 1):

        wo = wo.norm()
        wi = wi.norm()
        wh = (wo+wi).norm()

        rotation = matrix.Rotate(wh)
        wor = rotation.rotate(wo)
        wir = rotation.rotate(wi)
        wop = wor
        wip = wir
        wop.y = 0
        wip.y = 0

        wop = wop.norm()
        wip = wip.norm()
        thetaO = computeTheta(wop, 0)
        thetaI = computeTheta(wip, 0)

       
        #ZIPIN CONVENTION SAYS CHI PHI ON THE SAME SIDE IS NEGATIVE 
        ##place to change angle convention
        if thetaI * thetaO > 0:
            thetaI = math.fabs(thetaI) *  (-1)
        else:
            thetaI = math.fabs(thetaI)

 
        thetaO = math.fabs(thetaO)
        diff = thetaI - thetaO
        if math.fabs(diff) < 1e-5:
            thetaI = thetaO + 1e-5
            print (thetaO, thetaI) 
        

        #print(math.degrees(thetaO), math.degrees(thetaI))

        value = 0
        pdf = 0

        # for n == 1, we can skip far check if near check is successful
        for n in range(minBounce, maxBounce+1):
             
            ((nearHit, nearvalue, nearpdf), neartheta) = self.EvalNearHits(wor, wir, thetaO, thetaI, n, wop, wip)
            ((farHit, farvalue, farpdf), fartheta) = self.EvalFarHits(wor, wir, thetaO, thetaI, n, wop, wip )
            value += nearvalue + farvalue
            pdf += nearpdf + farpdf
            
        return (value, pdf)
                
            




        
