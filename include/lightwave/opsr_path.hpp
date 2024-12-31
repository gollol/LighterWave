#pragma once

#include <lightwave/opsr_data.hpp>

namespace lightwave {

    // Fast integer exponentiation
    // https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
    static int ipow(int base, int exp) {
        int result = 1;
        for (;;) {
            if (exp & 1)
                result *= base;
            exp >>= 1;
            if (!exp)
                break;
            base *= base;
        }
        return result;
    }

    // quantizationLevel is basicly whether its the 2d bin or the 3d or 4d or 5d bin (which implicitly is the bin size as well)
    static int roughnessToBin_weird(float roughness, int quantizationLevel) {
        return (int) std::min(
            std::floor(
                (std::pow(2.f, std::sqrt(roughness)) - 1) * (quantizationLevel + 1))
        , (float) quantizationLevel - 1);
    }

    // quantizationLevel is basicly whether its the 2d bin or the 3d or 4d or 5d bin (which implicitly is the bin size as well)
    static int roughnessToBin(float vertexRoughness, float roughnessRes) {
        float logRoughnessRes = roughnessRes + 1;
        return (int) std::min(
            std::floor(
                (std::pow(2.0f, std::sqrt(vertexRoughness)) - 1.0f) * logRoughnessRes)
        , (float) roughnessRes - 1);
    }


    static float getAttenuationFactor5d(float roughnessRes, int size, int startIndex, float virtualRoughness, 
    std::vector<float> &pathRoughness, AttenuationFactors attFactors) {

        int maxSubpathLength = 5;
        int subpathLength = min(size, maxSubpathLength);
        
        int queryIndex = 0;

        for (int i = startIndex; i < startIndex + subpathLength; i++) {
            if (i == startIndex && virtualRoughness != -1.0f) {
                // Use the virtual roughness instead of the start_idx for the
                // first vertex if a virtual roughness was accumulated
                queryIndex += roughnessToBin(virtualRoughness, roughnessRes);
            } else {
                queryIndex +=
                    ipow((int) roughnessRes, (i - startIndex)) *
                    roughnessToBin(pathRoughness[i], roughnessRes);
            }
        }

        float attenuationFactor = 0.0f;

        switch (subpathLength) {
            case 2:
                attenuationFactor =
                    attFactors.d2[queryIndex];
                break;
            case 3:
                attenuationFactor =
                    attFactors.d3[queryIndex];
                break;
            case 4:
                attenuationFactor =
                    attFactors.d4[queryIndex];
                break;
            case 5:
                attenuationFactor =
                    attFactors.d5[queryIndex];
                break;
            default:
                printf("Invalid subpath length %d", subpathLength);
        }

        return attenuationFactor;
    }

    static float accumulatedRoughnessFunction(std::vector<float> &pathRoughness) {
        float attenuationFactor = 0.5;
        float res = 0;

        res = 1 - pathRoughness.back();
        for (int i = 0; i < int(pathRoughness.size()) - 1; i++) {
           res *= 1 - attenuationFactor * pathRoughness[i];
        }
        res = 1 - res;
        return res;
    }

    static float getRougheningOPSR(float roughnessRes, std::vector<float> pathRoughness, AttenuationFactors attFactors) {
        // quick path for direct illumination
        if (pathRoughness.size() == 1)
            return pathRoughness.back();

        int processedIndex = 0;

        // Attenuation factor for the first subpath
        float attenuationFactor = getAttenuationFactor5d(roughnessRes, pathRoughness.size(), 0, -1.0f, pathRoughness, attFactors);

        // Keep the first vertex in the path separately to update it with the
        // virtual roughness for longer paths
        float virtualRoughness = pathRoughness[0];

        // Recursively estimate the attenuated roughening in the path
        while (true) {
            float accumulatedRoughness = 1.0f;

            if ((pathRoughness.size() - processedIndex) <= 5) // base case
            {
                // Accumulate the roughness with a potential virtual roughness
                accumulatedRoughness *=
                    (1.0f - attenuationFactor * virtualRoughness);
                for (int i = processedIndex + 1; i < int(pathRoughness.size()) - 1; i++) {
                    accumulatedRoughness *=
                        (1.0f - attenuationFactor * pathRoughness[i]);
                }
                accumulatedRoughness *= (1.0f - pathRoughness.back());

                return 1.0f - accumulatedRoughness;
            } else {
                // accumulate 5 vertex and update path roughness with updated
                // virtual roughness
                int lastIndex = processedIndex + 4;
                for (int i = processedIndex; i < lastIndex; i++) {
                    accumulatedRoughness *=
                        (1.0f - attenuationFactor * pathRoughness[i]);
                }
                accumulatedRoughness *=
                    (1.0f - pathRoughness[lastIndex]);

                // update processed index (5 - the virtual roughness updated)
                processedIndex += 4;

                // update the virtual roughness of the new subpath
                virtualRoughness = (1.0f - accumulatedRoughness);

                // query next attenuation factor
                attenuationFactor = getAttenuationFactor5d(
                    roughnessRes, pathRoughness.size() - processedIndex, lastIndex,
                    (1.0f - accumulatedRoughness), pathRoughness, attFactors);
            }
        }
    }



}