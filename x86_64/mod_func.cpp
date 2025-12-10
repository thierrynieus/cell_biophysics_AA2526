#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _Grc_sine_reg(void);
extern void _leak_reg(void);
extern void _potassium_reg(void);
extern void _ramp_reg(void);
extern void _sodium_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"Grc_sine.mod\"");
    fprintf(stderr, " \"leak.mod\"");
    fprintf(stderr, " \"potassium.mod\"");
    fprintf(stderr, " \"ramp.mod\"");
    fprintf(stderr, " \"sodium.mod\"");
    fprintf(stderr, "\n");
  }
  _Grc_sine_reg();
  _leak_reg();
  _potassium_reg();
  _ramp_reg();
  _sodium_reg();
}

#if defined(__cplusplus)
}
#endif
