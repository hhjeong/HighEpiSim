#ifndef __RNG__
#define __RNG__
#include <ctime>

// WELL512 random number generator

#define __W 32
#define __R 16
#define __P 0
#define __M1 13
#define __M2 9
#define __M3 5

#define __MAT0POS(t,v) (v^(v>>t))
#define __MAT0NEG(t,v) (v^(v<<(-(t))))
#define __MAT3NEG(t,v) (v<<(-(t)))
#define __MAT4NEG(t,b,v) (v ^ ((v<<(-(t))) & b))

#define __V0            __STATE[__state_i                   ]
#define __VM1           __STATE[(__state_i+__M1) & 0x0000000fU]
#define __VM2           __STATE[(__state_i+__M2) & 0x0000000fU]
#define __VM3           __STATE[(__state_i+__M3) & 0x0000000fU]
#define __VRm1          __STATE[(__state_i+15) & 0x0000000fU]
#define __VRm2          __STATE[(__state_i+14) & 0x0000000fU]
#define __newV0         __STATE[(__state_i+15) & 0x0000000fU]
#define __newV1         __STATE[__state_i                 ]
#define __newVRm1       __STATE[(__state_i+14) & 0x0000000fU]

#define __FACT 2.32830643653869628906e-10

static unsigned int __state_i = 0;
static unsigned int __STATE[__R]; 


static unsigned int __z0, __z1, __z2;

void InitWELLRNG512a () {
    int j;
    __state_i = 0;
    for (j = 0; j < __R; j++)
        __STATE[j] = (unsigned int)time(NULL);
}

double Random() {
    __z0    = __VRm1;
    __z1    = __MAT0NEG (-16,__V0)    ^ __MAT0NEG (-15, __VM1);
    __z2    = __MAT0POS (11, __VM2)  ;
    __newV1 = __z1                  ^ __z2; 
    __newV0 = __MAT0NEG (-2,__z0)     ^ __MAT0NEG(-18,__z1)    ^ __MAT3NEG(-28,__z2) ^ __MAT4NEG(-5,0xda442d24U,__newV1) ;
    __state_i = (__state_i + 15) & 0x0000000fU;
    return ((double) __STATE[__state_i]) * __FACT;
}

int RandomI( int mini, int maxi ) {
    // Output random integer in the interval min <= x <= max
    // Relative error on frequencies < 2^-32
    if (maxi <= mini) {
        if (maxi == mini) return mini; else return 0x80000000;
    }
    // Multiply interval with random and truncate
    int r = int((double)(unsigned int)(maxi - mini + 1) * Random() + mini); 
    if (r > maxi) r = maxi;
    return r;
}

#endif