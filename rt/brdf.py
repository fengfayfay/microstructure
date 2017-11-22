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


class Brdf:
    def __init__(self, alpha_x, alpha_y):
        self.microfacet = microfacet.Microfacet(alpha_x, alpha_y)
   
    def Pdf(self, wo, wh):
        pdf = self.microfacet.Pdf(wh)
        return pdf/(4.0* vec3.dot(wo, wh))


    #given wo, return sampled wi

    def MicrofacetValue(self, wo, wi, wh):
        value = self.microfacet.D(wh)
        cosThetaI = microfacet.CosTheta(wi)
        cosThetaO = microfacet.CosTheta(wo)
        value *= .25/(cosThetaI*cosThetaO)
        return value
        
    def Sample(self, wo, u):
        wo = wo.norm()
        (wh, cosTheta, phi) = self.microfacet.Sample_wh(wo, u)
        wi = ray.reflect(wh, -wo)
        value = self.MicrofacetValue(wo, wi, wh)
        pdf = self.microfacet.Pdf(wh) * .25 / vec3.dot(wo, wh)
        return (value, pdf, wi)
    

    def Eval(self, wo, wi):
        wo = wo.norm()
        wh = (wo + wi).norm()
        value = self.MicrofacetValue(wo, wi, wh)
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
        (wh, cosTheta, phi) = self.microfacet.Sample_wh(wo, u)
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
        pdf = self.Pdf(wo, wh) * prob

        thetaH = sample[0]
        #each sample contains the exit angle and the bounce count
        dprint(sample)
        sinThetaH = math.sin(thetaH)
        cosThetaH = math.cos(thetaH)
        wi = microfacet.SphericalDirection(sinThetaH, cosThetaH, phi)
        value = self.MicrofacetValue(wo, wi, wh) * prob
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
        microfacetD = self.MicrofacetValue(wo, wi, wh)
        microfacetPdf = self.Pdf(wi, wh)
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
            if nearHit:
                value += tmpvalue
                pdf += tmppdf
            else:
                (farHit, tmpvalue, tmppdf) = self.EvalFarHits(wo, wi, groovePhi, grooveChi, n)
                value += tmpvalue
                pdf += tmppdf
        return (value, pdf)
                
            




        
