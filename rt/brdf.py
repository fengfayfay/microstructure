import microfacet
import zipinPaper
import math
import vec3
import ray

#given wo, return sampled wi direction

class Brdf:
    def __init__(self, alpha_x, alpha_y):
        self.microfacet = microfacet.Microfacet(alpha_x, alpha_y)
   
    def Pdf(self, wo, wh):
        pdf = self.microfacet.Pdf(wh)
        return pdf/(4.0* vec3.dot(wo, wh))


    #given wo, return sampled wi
    def Sample(self, wo, u):
        wo = wo.norm()
        (wh, cosTheta, phi) = self.microfacet.Sample_wh(wo, u)
        wi = ray.reflect(wh, -wo)
        value = self.microfacet.D(wh)
        cosThetaI = microfacet.CosTheta(wi)
        cosThetaO = microfacet.CosTheta(wo)
        value *= .25/(cosThetaI*cosThetaO)
        pdf = self.microfacet.Pdf(wh) * .25 / vec3.dot(wo, wh)
        return (value, pdf, wi)
    

    def Eval(self, wo, wi):
        wo = wo.norm()
        wh = (wo + wi).norm()
        value = self.microfacet.D(wh)
        cosThetaI = microfacet.CosTheta(wi)
        cosThetaO = microfacet.CosTheta(wo)
        value *= .25/(cosThetaI*cosThetaO)
        pdf = self.Pdf(wo, wh)
        return (value, pdf)


        
