// For use in counter_METTriggers.py

//float normalizeby_nGoodPV_contents[40] = {0.3946, 0.0000, 2.2613, 1.3845, 1.7964, 1.2426, 1.4753, 1.3611, 1.5793, 1.3593, 1.4691, 1.2445, 1.0870, 1.0926, 0.9985, 1.0597, 0.8314, 0.8618, 0.9619, 0.8383, 0.8900, 0.8989, 0.7458, 0.6986, 0.9362, 0.7171, 1.2954, 1.0749, 2.1177, 1.0999, 0.9589, 1.0307, 3.0766, 0.6461, 0.7384, 0.3692, 0.8614, 0.7384, 0.2307, 1.4460};
//float normalizeby_nGoodPV_contents[40] = {0.3946, 2.2613, 2.2613, 1.3845, 1.7964, 1.2426, 1.4753, 1.3611, 1.5793, 1.3593, 1.4691, 1.2445, 1.0870, 1.0926, 0.9985, 1.0597, 0.8314, 0.8618, 0.9619, 0.8383, 0.8900, 0.8989, 0.7458, 0.6986, 0.9362, 0.7171, 1.1946, 1.1946, 1.5361, 1.5361, 0.9876, 0.9876, 0.9037, 0.9037, 0.7015, 0.7015};  // merge some bins

//float normalizeby_noiseType_contents[5] = {1.0428, 1.0009, 0.9747, 1.3487, 0.2501};
//float normalizeby_noiseType_contents[5] = {1.0478, 0.9991, 0.9940, 1.3229, 0.2553};

//float normalizeby_noiseType2_contents[3] = {1.0428, 1.0009, 0.8763};
//float normalizeby_noiseType2_contents[3] = {1.0478, 0.9991, 0.8888};

//float normalizeby_nGoodPV_contents[40] = {0.3338, 0.0000, 0.0000, 0.0246, 0.0364, 0.0502, 0.0922, 0.1161, 0.2171, 0.2733, 0.4338, 0.5130, 0.6587, 0.8274, 0.9798, 1.2872, 1.1953, 1.4292, 1.7392, 1.6665, 1.8842, 2.0042, 1.7360, 1.7080, 2.3386, 1.8115, 3.2579, 2.7722, 5.5341, 2.8798, 2.1785, 2.7055, 8.1986, 1.6944, 1.9677, 0.9838, 2.2956, 1.9677, 0.6149, 3.8533};
float normalizeby_nGoodPV_contents[40] = {0.3338, 0.0189, 0.0189, 0.0246, 0.0364, 0.0502, 0.0922, 0.1161, 0.2171, 0.2733, 0.4338, 0.5130, 0.6587, 0.8274, 0.9798, 1.2872, 1.1953, 1.4292, 1.7392, 1.6665, 1.8842, 2.0042, 1.7360, 1.7080, 2.0750, 2.0750, 3.0358, 3.0358, 4.0173, 4.0173, 2.3701, 2.3701, 3.3204, 3.3204, 1.5634, 1.5634, 1.5634, 1.5634, 1.6944, 1.6944};

float normalizeby_noiseType_contents[5] = {0.9737, 1.0154, 0.9307, 1.2614, 0.2980};

float normalizeby_noiseType2_contents[3] = {0.9737, 1.0154, 0.8558};


inline float normalizeby_nGoodPV(float npv) {
    if (npv<=0.) return normalizeby_nGoodPV_contents[0];
    if (npv>40.) return normalizeby_nGoodPV_contents[39];
    int bin = npv;
    return normalizeby_nGoodPV_contents[bin];
}

inline float normalizeby_noiseType(float f0, float f1, float f2, float f3) {
    if (int(f0)==0)
        return normalizeby_noiseType_contents[1];
    else if (int(f1)==0)
        return normalizeby_noiseType_contents[2];
    else if (int(f2)==0)
        return normalizeby_noiseType_contents[3];
    else if (int(f3)==0)
        return normalizeby_noiseType_contents[4];
    return normalizeby_noiseType_contents[0];
}

inline float normalizeby_noiseType2(float f0, float f1) {
    if (int(f0)==0)
        return normalizeby_noiseType2_contents[1];
    else if (int(f1)==0)
        return normalizeby_noiseType2_contents[2];
    return normalizeby_noiseType2_contents[0];
}
