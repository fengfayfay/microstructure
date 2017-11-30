import microfacet
import zipinPaper
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
        self.microfacet = microfacet.Microfacet(alpha_x, alpha_y)
   
    def Pdf(self, wo, wh):
        pdf = self.microfacet.Pdf(wh)
        return pdf/(4.0* vec3.dot(wo, wh))


    #given wo, return sampled wi

    def MicrofacetValue(self, wo, wi, wh, withoutG = False):
        cosThetaI = microfacet.CosTheta(wi)
        cosThetaO = microfacet.CosTheta(wo)
        if cosThetaI == 0 or cosThetaO == 0:
            return 0
        if wh.x==0 and wh.y ==0 and wh.z == 0:
            return 0
        value = self.microfacet.D(wh)
        #value *= .25/(cosThetaI*cosThetaO)
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
        #if wi.z * wo.z <= 0:
        #    return (0, 0, None)
        value = self.MicrofacetValue(wo, wi, wh, withoutG)
        pdf = self.microfacet.Pdf(wh) * .25 / vec3.dot(wo, wh)
        return (value, pdf, wi)
    

    def Eval(self, wo, wi, withoutG = False):
        wo = wo.norm()
        wh = (wo + wi).norm()
        value = self.MicrofacetValue(wo, wi, wh, withoutG)
        pdf = self.Pdf(wo, wh)
        return (value, pdf)



ANGLEERROR = .00001
MAXBOUNCE = 20

class ZipinBrdf(Brdf):

    #TODO:  figure out the correct normalization factors
    #       and angle orientations for grooveTheta
    #       groovePhi
    #       and final grooveChi
    def Sample(self, wo, u):
        wo = wo.norm()
        if wo.z == 0:
            return (0, 0, None, None)
        (wh, cosTheta, phi) = self.microfacet.Sample_wh(wo, u)
        if wh.x == 0 and wh.y == 0 and wh.z == 0:
            return (0, 0, None, None)
        grooveAlpha = math.acos(cosTheta)
        grooveTheta = math.pi * .5 - grooveAlpha 
        cosThetaO = microfacet.CosTheta(wo)
        groovePhi = math.acos(cosThetaO)

        
        hits = zipinPaper.zipinPaper(grooveTheta, groovePhi, False)
        dprint(hits)
        #dict = collections.defaultdict(float)
        #for hit in hits:
        #    dict[hit[0]] = hit[1]
        (sample, prob) = weight.weightedRandomChoice(hits)

        thetaH = sample[0]
        #each sample contains the exit angle and the bounce count
        dprint(sample)
        sinThetaH = math.sin(thetaH)
        cosThetaH = math.cos(thetaH)
        wi = microfacet.SphericalDirection(sinThetaH, cosThetaH, phi)

        scaleD = .25 * vec3.dot(wi, wh)
        microfacetD = self.microfacet.D(wh) * scaleD
        scalePdf = .25 * vec3.dot(wo, wh)
        microfacetPdf = self.microfacet.Pdf(wh) * scalePdf
        value = microfacetD * prob
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
        wh = microfacet.SphericalDirection(math.sin(grooveAlpha), math.cos(grooveAlpha), random.uniform(0, 1) * math.pi * 2.0)
        scaleD = .25 * vec3.dot(wi, wh)
        microfacetD = self.microfacet.D(wh) * scaleD
        scalePdf = .25 * vec3.dot(wo, wh)
        microfacetPdf = self.microfacet.Pdf(wh) * scalePdf
        #dprint(wh, microfacetD, microfacetPdf, side)
        hits = zipinPaper.zipinPaper(grooveTheta, groovePhi, False)
       
        value = 0;
        pdf = 0;
        for hit in hits:
            if hit[3] != side:
                continue
            hitChi = hit[0][0]
            diff = math.fabs(math.fabs(hitChi) - math.fabs(grooveChi))
            #diff = math.fabs(hitChi - grooveChi)
            if diff < ANGLEERROR and hit[1] > 0:
                dprint(hit)
                value += microfacetD * hit[1]
                pdf += microfacetPdf * hit[1]
                return (True, value, pdf)
        return (False, 0, 0)
        

    def Eval(self, wo, wi):
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

        for n in range(1, MAXBOUNCE):
            (nearHit, tmpvalue, tmppdf) = self.EvalNearHits(wo, wi, groovePhi, grooveChi, n)
            value += tmpvalue
            pdf += tmppdf
            (farHit, tmpvalue, tmppdf) = self.EvalFarHits(wo, wi, groovePhi, grooveChi, n)
            value += tmpvalue
            pdf += tmppdf
        return (value, pdf)
                
            




        
