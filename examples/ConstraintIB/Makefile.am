## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = eel2d eel3d externalFlowBluffBody2d externalFlowBluffBody2d_stabilized

## Standard make targets.
examples:
	@(cd eel2d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
	@(cd eel3d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
  @(cd externalFlowBluffBody2d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
  @(cd externalFlowBluffBody2d_stabilized && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
