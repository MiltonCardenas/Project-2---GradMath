# -*- coding: utf-8 -*-
"""
@author: milton

AUXILIAR CODES
"""
import random

def vary_vector(vector, relative_difference):
    varied_vector = []
    for value in vector:
        variation = value * relative_difference
        random_variation = random.uniform(-variation, variation)
        varied_vector.append(value + random_variation)
    return varied_vector