ch.pdf: ch.tex cahn-hilliard.bib
	cp ch.tex tmp.tex
	latex tmp && bibtex tmp && latex tmp && latex tmp
	dvipdf tmp.dvi ch.pdf
	dvips tmp.dvi -o ch.ps

clean:
	rm -f tmp.* *.ps *.pdf
