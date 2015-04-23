BIB = bibtex
LTX = pdflatex
TOPDF = dvipdf
RUNPDF = xdg-open
MAIN = meta_strategy_model
SI = meta_strategy_model-si

all: main
main: main 
si: si

quick: ${MAIN}.aux ${MAIN}.bbl
	pdflatex ${MAIN}.tex

clean:
	rm -f ${MAIN}.aux 
	rm -f ${MAIN}.log 
	rm -f ${MAIN}.blg 
	rm -f ${MAIN}.bbl 
	rm -f ${MAIN}.dvi 
	rm -f ${MAIN}.out 
	rm -f ${MAIN}.pdf
        
	rm -f ${SI}.aux 
	rm -f ${SI}.log 
	rm -f ${SI}.blg 
	rm -f ${SI}.bbl 
	rm -f ${SI}.dvi 
	rm -f ${SI}.out 
	rm -f ${SI}.pdf

main: meta_strategy_model.tex 
	${LTX} ${MAIN}.tex
	${BIB} ${MAIN}
	${LTX} ${MAIN}.tex
	${LTX} ${MAIN}.tex
	${RUNPDF} ${MAIN}.pdf
	#${TOPDF} ${MAIN}.dvi


si: meta_strategy_model-si.tex 
	${LTX} ${SI}.tex
	${BIB} ${SI}
	${LTX} ${SI}.tex
	${LTX} ${SI}.tex
	#${TOPDF} ${SI}.dvi

