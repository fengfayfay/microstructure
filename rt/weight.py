import math
import os.path
import random

# Function: Pdf
# -------------------------
# Returns the probability density of a Gaussian distribution with
# the specified mean and std, evaluated at the specified value.
def pdf(mean, std, value):
    u = float(value - mean) / abs(std)
    y = (1.0 / (math.sqrt(2 * math.pi) * abs(std))) * math.exp(-u * u / 2.0)
    return y

# Function: Weighted Random Choice
# --------------------------------
# Given a dictionary of the form element -> weight, selects an element
# randomly based on distribution proportional to the weights. Weights can sum
# up to be more than 1. 
def weightedRandomChoice(weightDict, maxBounce = 1):
    weights = []
    elems = []
    for elem in weightDict:
#        if elem[2] > maxBounce:
#            continue
        weights.append(elem[1])
        elems.append(elem)
    #print(weights)
    total = sum(weights)
    key = random.uniform(0, total)
    runningTotal = 0.0
    chosenIndex = None
    for i in range(len(weights)):
        weight = weights[i]
        runningTotal += weight
        prob = weight/total
        #print(prob)
        if runningTotal > key:
            chosenIndex = i
            return (elems[chosenIndex], prob)
    return (None, 0)

