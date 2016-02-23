docfiles = manual.tex
if HAVE_PDFLATEX
doc_DATA = manual.pdf
MANNAME = manual
MANTEXSRC = $(MANNAME).tex coal_manual.tex
MANAUX = $(MANNAME).aux
MANPDF = $(MANNAME).pdf
DOCDIR=$(abs_top_builddir)/doc
CLEANFILES = $(MANPDF) $(MANNAME).log $(MANNAME).idx $(MANNAME).out $(MANNAME).toc $(MANAUX)

$(MANPDF):  $(DOCDIR)/$(MANTEXSRC)
	$(PDFLATEX) $(MANTEXSRC)
#	$(BIBTEX) $(MANAUX)
#	$(PDFLATEX) $(MANTEXSRC)
	$(PDFLATEX) $(MANTEXSRC)
dist_doc_DATA = $(docfiles)
endif

EXTRA_DIST = manual.pdf

#SUBDIRS =
