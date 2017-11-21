import microfacet
import zipinPaper
import math
import vec3
import ray
import weight
import collections
import random

#given wo, return sampled wi direction

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
        print(hits)
        #dict = collections.defaultdict(float)
        #for hit in hits:
        #    dict[hit[0]] = hit[1]
        (sample, prob) = weight.weightedRandomChoice(hits)
        pdf = self.microfacet.Pdf(wh) * prob

        thetaH = sample[0]
        #each sample contains the exit angle and the bounce count
        print(sample)
        sinThetaH = math.sin(thetaH)
        cosThetaH = math.cos(thetaH)
        wi = microfacet.SphericalDirection(sinThetaH, cosThetaH, phi)
        #value = self.MicrofacetValue(wo, wi, wh)
        value = self.microfacet.D(wh) * prob
        return (value, pdf, wi, wh)

    def Eval(self, wo, wi):
        wo = wo.norm()
        wi = wi.norm()
        cosThetaO = microfacet.CosTheta(wo)
        cosThetaI = microfacet.CosTheta(wi)
        groovePhi = math.acos(cosThetaO)
        grooveChi = math.acos(cosThetaI)  
        print("phi: " + repr(groovePhi))
        print("chi: "+ repr(grooveChi))

        value = 0
        pdf = 0

        for n in range(1, MAXBOUNCE):
            grooveTheta = math.fabs(math.pi - groovePhi - grooveChi) * .5 / n
            grooveAlpha = math.pi * .5 - grooveTheta
            wh = microfacet.SphericalDirection(math.sin(grooveAlpha), math.cos(grooveAlpha), random.uniform(0, 1) * math.pi * 2.0)
            microfacetD = self.microfacet.D(wh)
            microfacetPdf = self.microfacet.Pdf(wh)
            print(wh)
            print(microfacetD)
            print(microfacetPdf)
            
            hits = zipinPaper.zipinPaper(grooveTheta, groovePhi, False)
            #print(hits)
            #hit[0] is sample, hit[1] is weight
            for hit in hits:
                hitChi = hit[0][0]
                diff = math.fabs(math.fabs(hitChi) - math.fabs(grooveChi))
                if diff < ANGLEERROR and hit[1] > 0:
                    print(hit)
                    value += microfacetD * hit[1]
                    pdf += microfacetPdf * hit[1]
                    break
        return (value, pdf)
                
            




        
