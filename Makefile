.PHONY: doc



doc:
	sphinx-apidoc . -o doc/source/podi/


all: doc html latex

html: doc
	$(MAKE) html -C doc/

latex: doc
	$(MAKE) latex -C doc/
	cd doc/build/latex; pwd; pdflatex -interaction=nonstopmode QuickReduce.tex
	cp doc/build/latex/QuickReduce.pdf doc/
