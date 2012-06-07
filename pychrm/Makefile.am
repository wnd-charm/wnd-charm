
BUILT_SOURCES = $(srcdir)/swig/pychrm_wrap.cxx
SWIG_SOURCES = $(srcdir)/swig/pychrm.i \
	$(srcdir)/swig/cmatrix.i \
	$(srcdir)/swig/wndchrm_error.i \
	$(srcdir)/swig/FeatureAlgorithms.i \
	$(srcdir)/swig/FeatureTransforms.i

pkgpython_PYTHON = pychrm.py TrainingSet.py FeatureNameMap.py FeatureRegistration.py __init__.py
pkgpyexec_LTLIBRARIES = _pychrm.la
_pychrm_la_SOURCES = $(srcdir)/swig/pychrm_wrap.cxx $(SWIG_SOURCES)
_pychrm_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_pychrm_la_LDFLAGS = -module -bundle 
_pychrm_la_LIBADD = ../src/libchrm.la

$(srcdir)/swig/pychrm_wrap.cxx : $(SWIG_SOURCES)
	 $(SWIG) $(SWIG_PYTHON_OPT) -I$(top_srcdir)/src -outdir $(srcdir) -o $@ $<


