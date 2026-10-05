// Internal unit system for the whole integrator: length in kpc, time in Myr,
// velocity in kpc/Myr, mass in solar masses. Fully self-consistent (no hidden
// factors inside force/integrator code) -- convert at the edges only, with the
// helpers below, when data comes in as km/s (the usual astronomical unit).
#ifndef NOSTOS_UNITS_H
#define NOSTOS_UNITS_H

// 1 kpc / 1 Myr, expressed in km/s. kpc->km and Myr->s are both exact SI
// definitions (Julian year = 365.25 days), so this is exact to the precision
// shown, not a fitted/rounded value.
//   1 kpc = 3.0856775814913673e16 km
//   1 Myr = 3.15576e13 s   (365.25 * 86400 * 1e6)
//   (km/s) per (kpc/Myr) = 3.0856775814913673e16 / 3.15576e13
#define NOSTOS_KMS_PER_KPCMYR 977.7922216807892

static inline double nostos_kms_to_kpcmyr(double v_kms) {
    return v_kms / NOSTOS_KMS_PER_KPCMYR;
}
static inline double nostos_kpcmyr_to_kms(double v_kpcmyr) {
    return v_kpcmyr * NOSTOS_KMS_PER_KPCMYR;
}

// Newton's constant, G = 4.30091e-6 kpc*(km/s)^2 / Msun (standard value quoted
// e.g. by galpy/astropy; derivable from G = 6.6743e-11 m^3 kg^-1 s^-2 plus
// Msun = 1.98892e30 kg and 1 kpc = 3.0856775814913673e19 m).
// Converted here into the internal unit system: kpc^3 Msun^-1 Myr^-2.
#define NOSTOS_G_KPC_KMS2_PER_MSUN 4.30091e-6
#define NOSTOS_G (NOSTOS_G_KPC_KMS2_PER_MSUN / (NOSTOS_KMS_PER_KPCMYR * NOSTOS_KMS_PER_KPCMYR))

#endif
