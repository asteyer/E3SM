NAME=FirstnameLastname

all: ${NAME}.pdf

${NAME}.pdf: ${NAME}.tex ${NAME}.bib
	pdflatex ${NAME}
	bibtex ${NAME}
	pdflatex ${NAME}
	pdflatex ${NAME}

clean:
	rm ${NAME}.pdf ${NAME}.bbl ${NAME}.blg ${NAME}.log ${NAME}.aux
