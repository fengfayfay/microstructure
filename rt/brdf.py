import microfacet
import zipinPaper
import zipin
import math
import vec3
import ray
import weight
import collections
import random

#given wo, return sampled wi direction
debug = 0

def dprint(s):
    if debug == 1:
        print(s)

def reflect(wo, n):
    return -wo + n * vec3.dot(wo, n) * 2.0


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

    def MicrofacetValue(self, wo, wi, wh, withoutG = False):
        cosThetaI = microfacet.CosTheta(wi)
        cosThetaO = microfacet.CosTheta(wo)
        if cosThetaI == 0 or cosThetaO == 0:
            return 0
        if wh.x==0 and wh.y ==0 and wh.z == 0:
            return 0
        value = self.microfacet.D(wh)
        value *= .25/(cosThetaO)
        if withoutG:
            return value
        else:
            return  value * self.microfacet.G(wo, wi, wh)
        
    def Sample(self, wo, u, withoutG=False):
        wo = wo.norm()
        if wo.z == 0:
            return (0, 0, None)
        (wh, cosTheta, phi) = self.microfacet.Sample_wh(wo, u)
        if wh.x == 0 and wh.y == 0 and wh.z == 0:
            return (0, 0, None)
        wi = reflect(wo, wh)
        value = self.MicrofacetValue(wo, wi, wh, withoutG)
        pdf = self.Pdf(wo, wh) 
        return (value, pdf, wi)
    

    def Eval(self, wo, wi, withoutG = False):
        wo = wo.norm()
        wh = (wo + wi).norm()
        value = self.MicrofacetValue(wo, wi, wh, withoutG)
        pdf = self.Pdf(wo, wh)
        return (value, pdf)



ANGLEERROR = .00001
MAXBOUNCE = 2 
MICROFACETTHRESH = .000001

class ZipinBrdf(Brdf):
    def __init__(self, alpha_x, alpha_y, zipinVersion):
        self.microfacet = microfacet.Beckmann(alpha_x, alpha_y)
        print('zipinVersion: ', zipinVersion)
        if zipinVersion == 'Feng':
            self.zipin = zipinPaper.zipin
        else:
            self.zipin = zipin.zipin

    #TODO:  figure out the correct normalization factors
    #       and angle orientations for grooveTheta
    #       groovePhi
    #       and final grooveChi
    def Sample(self, wo, u):
        wo = wo.norm()
        if wo.z == 0:
            return (0, 0, None, None)
        (wh, cosThetaH, phiH) = self.microfacet.Sample_wh(wo, u)
        if wh.x == 0 and wh.y == 0 and wh.z == 0:
            return (0, 0, None, None)
        grooveAlpha = math.acos(cosThetaH)
        grooveTheta = math.pi * .5 - grooveAlpha 
        cosThetaO = microfacet.CosTheta(wo)
        groovePhi = math.acos(cosThetaO)

        
        hits = self.zipin(math.degrees(grooveTheta), math.degrees(groovePhi))
        (sample, prob) = weight.weightedRandomChoice(hits)

        thetaH = math.radians(sample)
        dprint(sample)
        sinThetaH = math.sin(thetaH)
        cosThetaH = math.cos(thetaH)
        wi = microfacet.SphericalDirection(sinThetaH, cosThetaH, phiH)

        microfacetValue = self.MicrofacetValue(wo, wi, wh, True)
        microfacetPdf = self.Pdf(wo, wh)
        value = microfacetValue * prob
        pdf = microfacetPdf * prob

        return (value, pdf, wi, wh)

    def EvalNearHits(self, wo, wi, phi, chi, n):
        sign = -1 if n % 2 == 0 else 1
        chi *= sign
        theta = (math.pi + phi - chi) *.5/n
        return self.EvalHits(wo, wi, chi, phi, theta, n, 'right')

    def EvalFarHits(self, wo, wi, phi, chi, n):
        sign = -1 if n % 2 else 1
        chi *= sign
        theta = (math.pi -phi -chi) * .5/n
        return self.EvalHits(wo, wi, chi, phi, theta, n, 'left')

    def EvalHits(self, wo, wi, grooveChi, groovePhi, grooveTheta, n, side):
         
        dprint("theta: "+ repr(grooveTheta))
        grooveAlpha = math.pi * .5 - grooveTheta
        if side == 'left':
            wh = microfacet.SphericalDirection(-math.sin(grooveAlpha), math.cos(grooveAlpha), random.uniform(0, 1) * math.pi * 2.0)
        else:
            wh = microfacet.SphericalDirection(math.sin(grooveAlpha), math.cos(grooveAlpha), random.uniform(0, 1) * math.pi * 2.0)

        microfacetValue = self.MicrofacetValue(wo, wi, wh, True)
        microfacetPdf = self.Pdf(wo, wh)
        if microfacetValue < MICROFACETTHRESH:
            return (False, 0, 0)
        hits = self.zipin(math.degrees(grooveTheta), math.degrees(groovePhi))
       
        value = 0;
        pdf = 0;
        for hit in hits:
            if hit[3] != side:
                continue
            hitChi = math.radians(hit[0])
            #diff = math.fabs(math.fabs(hitChi) - math.fabs(grooveChi))
            diff = math.fabs(hitChi - grooveChi)
            if diff < ANGLEERROR and hit[1] > 0:
                #print(groovePhi, grooveTheta, wh, microfacetD, microfacetPdf, hit[1], side)
                #print(hit)
                value += microfacetValue * hit[1]
                pdf += microfacetPdf * hit[1]
        return (value>0, value, pdf)
        

    def Eval(self, wo, wi, maxBounce = MAXBOUNCE):
        wo = wo.norm()
        wi = wi.norm()
        cosThetaO = microfacet.CosTheta(wo)
        cosThetaI = microfacet.CosTheta(wi)
        groovePhi = math.acos(cosThetaO)
        grooveChi = math.acos(cosThetaI)  
        dprint("phi: " + repr(groovePhi))
        dprint("chi: "+ repr(grooveChi))

        value = 0
        pdf = 0

        for n in range(maxBounce):
            (nearHit, tmpvalue, tmppdf) = self.EvalNearHits(wo, wi, groovePhi, grooveChi, n+1)
            value += tmpvalue
            pdf += tmppdf
            #if nearHit:
            #    continue
            (farHit, tmpvalue, tmppdf) = self.EvalFarHits(wo, wi, groovePhi, grooveChi, n+1)
            value += tmpvalue
            pdf += tmppdf
        return (value, pdf)
                
            




        
