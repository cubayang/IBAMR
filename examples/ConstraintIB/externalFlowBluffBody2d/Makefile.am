## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = IBRigidBodyKinematics.cpp main.cpp

EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += externalFlowBluffBody2d
endif

externalFlowBluffBody2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
externalFlowBluffBody2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
externalFlowBluffBody2d_SOURCES = $(SOURCES)

examples: $(EXTRA_PROGRAMS)
