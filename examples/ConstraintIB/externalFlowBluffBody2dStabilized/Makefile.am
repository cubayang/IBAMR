## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = IBRigidBodyKinematics.cpp main.cpp
EXTRA_DIST = input2d

EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += externalFlowBluffBody2dStabilized
endif

externalFlowBluffBody2dStabilized_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
externalFlowBluffBody2dStabilized_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
externalFlowBluffBody2dStabilized_SOURCES = $(SOURCES)

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		cp -f $(srcdir)/input2d $(PWD) ; \
		cp -f $(srcdir)/bluffBody2d_ds004.vertex $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input2d ; \
	  rm -f $(builddir)/bluffBody2d_ds004.vertex ; \
	fi ;